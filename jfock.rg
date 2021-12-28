import "regent"

require "helper"
require "md.jfock.generate_jfock_integral"

local c = regentlib.c

-----------------------------------------------------------------------------------------
-- Terra function that:
--  (1) generates Regent integral tasks of all angular momentum types that are required
--      to compute the J matrix up to a user-specified max angular momentum 
--  (2) executes the integral tasks at runtime using input data (regions of data for bras 
--      and kets and parameters for Schwartz screening threshold and desired parallelism)
-----------------------------------------------------------------------------------------
function jfock(r_jbras_list, r_jkets_list, r_gamma_table, threshold, parallelism)
  local statements = terralib.newlist()
  -- Generate all integral tasks up to and including largest desired angular momentum
  for L1 = 0, getCompiledMaxMomentum() do -- inclusive
    -- We generate one extra kernel so that we can compute the gradient.
    for L2 = L1, getCompiledMaxMomentum() + 1 do -- inclusive
      for L3 = 0, getCompiledMaxMomentum() do -- inclusive
        for L4 = L3, getCompiledMaxMomentum() do -- inclusive
          local r_jbras = r_jbras_list[L1][L2]
          local r_jkets = r_jkets_list[L3][L4]
          local L12, L34 = L1 + L2, L3 + L4
          -- Generate task (compile time)
          local jfock_integral = generateTaskMcMurchieJFockIntegral(L12, L34)
          -- only execute tasks if regions are not empty
          if r_jbras ~= nil and r_jkets ~= nil then
            statements:insert(rquote
              -- partition bra regions into ‘parallelism’ number of equally-sized regions
	            var jbra_coloring = ispace(int1d, parallelism)
	            var p_jbras = partition(equal, r_jbras, jbra_coloring)
              -- execute tasks (runtime)
              __demand(__index_launch)
              for jbra_color in jbra_coloring do
                jfock_integral(p_jbras[jbra_color], r_jkets, r_gamma_table, threshold)
              end
            end)
          end
        end
      end
    end
  end
  return statements
end
