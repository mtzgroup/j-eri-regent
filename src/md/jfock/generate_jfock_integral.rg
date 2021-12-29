import "regent"

require "fields"
require "helper"
require "md.generate_R_table"
require "md.jfock.generate_kernel_statements"

local rsqrt = regentlib.rsqrt(double)

-----------------------------------------------------------------------------------------
-- Terra function that, given angular momenta for the bra (L12) and ket (L34), returns 
--   the task ‘jfock_integral’ which computes J matrix elements using the McMurchie-   
--   Davidson algorithm                                                                
-----------------------------------------------------------------------------------------
local _jfock_integral_cache = {}
function generateTaskMcMurchieJFockIntegral(L12, L34)
  local H12, H34 = tetrahedral_number(L12 + 1), tetrahedral_number(L34 + 1)
  local L_string = LPairToStr[L12]..LPairToStr[L34]
  if _jfock_integral_cache[L_string] ~= nil then
    return _jfock_integral_cache[L_string]
  end
  -- Create an empty table of Regent variables ‘R’ to hold auxiliary values RNLMj.
  local R = {}
  for N = 0, L12+L34 do -- inclusive
    R[N] = {}
    for L = 0, L12+L34-N do -- inclusive
      R[N][L] = {}
      for M = 0, L12+L34-N-L do -- inclusive
        R[N][L][M] = {}
        for j = 0, L12+L34-N-L-M do -- inclusive
          -- Create a new Regent variable named RNLMj (where N,L,M,j are loop variables,
          -- e.g. ‘R1001’) and store it in the corresponding location in table ‘R’
          R[N][L][M][j] = regentlib.newsymbol(double, "R"..N..L..M..j)
        end
      end
    end
  end

  ----------------------------------------------------------------------------------
  -- Regent task that computes J matrix elements in Hermite basis: 
  --     J_[P] = \sum_[Q] [P|Q] * density_[Q]  (eq. 15 in main text)
  --     screened below the param ‘threshold’  (eq. 33 in main text)
  --
  -- For example, for an SSSP task (L12 = 0 and L34 = 1) with N bras (indexed by n) 
  -- and M kets (indexed by m):
  --   for all n, m 
  --       r_jkets[m].density contains intput density matrix values and is size 
  --            H34 = tetrahedral_number(2) = 4 
  --       r_jbras[n].output  contains output Coulomb matrix values and is size 
  --            H12 = tetrahedral_number(1) = 1 
  --    \sum_[Q] is the loop over the M kets and is performed for each of the N bras
  ----------------------------------------------------------------------------------
  local -- jfock_integal task is a local variable within the Terra function
  __demand(__leaf) -- tell compiler this task will call no other tasks (optimize)
  __demand(__cuda) -- tell compiler to generate GPU code and execute on a GPU
  task jfock_integral(r_jbras       : region(ispace(int1d), getJBra(L12)), -- bra basis and output region
                      r_jkets       : region(ispace(int1d), getJKet(L34)), -- ket basis and density region
                      r_gamma_table : region(ispace(int2d), double[5]),    -- tablulated data for Boys fn. interp.
                      threshold     : float)                               -- 4-index integral screening thresh.
  where -- declare region priviledges (used by runtime for task dependency analysis)
    reads(r_jbras.{location, eta, C, bound}, r_jkets, r_gamma_table),
    reduces +(r_jbras.output)
  do
    var jket_idx_bounds_lo : int = r_jkets.ispace.bounds.lo
    var jket_idx_bounds_hi : int = r_jkets.ispace.bounds.hi
    for jbra_idx in r_jbras.ispace do
      var jbra_location : Point = r_jbras[jbra_idx].location -- Point is a struct of (x,y,z)
      var jbra_eta : double = r_jbras[jbra_idx].eta
      var jbra_C : double = r_jbras[jbra_idx].C
      var jbra_bound : float = r_jbras[jbra_idx].bound
      var accumulator : double[H12] -- an array to hold output values 
      for i = 0, H12 do -- exclusive
        accumulator[i] = 0.0
      end

      for jket_idx = jket_idx_bounds_lo, jket_idx_bounds_hi + 1 do -- exclusive
        var jket = r_jkets[jket_idx]

        var bound : float = jbra_bound * jket.bound
        if bound <= threshold then break end -- Schwartz screening

        var a : double = jbra_location.x - jket.location.x
        var b : double = jbra_location.y - jket.location.y
        var c : double = jbra_location.z - jket.location.z

        var alpha : double = jbra_eta * jket.eta * (1.0 / (jbra_eta + jket.eta))
        var lambda : double = jbra_C * jket.C * rsqrt(jbra_eta + jket.eta)
        var t : double = alpha * (a*a + b*b + c*c);
        [generateStatementsComputeRTable(R, L12+L34+1, t, alpha, lambda,
                                         a, b, c, r_gamma_table)];

        [generateJFockKernelStatements(R, L12, L34, rexpr jket.density end,
                                       accumulator)]
      end
      r_jbras[jbra_idx].output += accumulator
    end
  end -- end jfock_integral task

  jfock_integral:set_name("JFockMcMurchie"..L_string) -- give this task a unique name so it can be called
  _jfock_integral_cache[L_string] = jfock_integral
  return jfock_integral
end
