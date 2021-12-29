import "regent"

require "helper"

----------------------------------------------------------------------------------------
-- Terra function that returns a list of Regent statements that combine elements of 
-- the density matrix (in the Hermite basis Q) with previously computed RNLM0 auxiliary 
-- functions (lambda already incorporated) to compute J_[P] (elements of J in the 
-- Hermite basis)
--
-- J_[P] = \sum_[Q] [P|Q] * density_[Q] 
-- [P|Q] = lambda * (-1)^(Nq+Lq+Mq) * R_(Np+Nq, Lp+Lq, Mp+Mq, 0)
-- (See eq. 13 and eq. 15 in main text)
--
-- For example, for an SSSP task (L12 = 0, H12 = 1 and L34 = 1, H34 = 4), this function
-- returns the following statement list:
--    results[0] : double = 0.0
--    results[0] += density[0] * R[0][0][0][0]
--    results[0] -= density[1] * R[1][0][0][0]
--    results[0] -= density[2] * R[0][1][0][0]
--    results[0] -= density[3] * R[0][0][1][0]
--    accumulator[0] += results[0]
----------------------------------------------------------------------------------------
function generateJFockKernelStatements(R, L12, L34, density, accumulator)
  local H12, H34 = tetrahedral_number(L12 + 1), tetrahedral_number(L34 + 1)
  -- create a Terra list onto which Regent statements will be appended
  local statements = terralib.newlist()
  local results = {}
  for i = 0, H12-1 do -- inclusive
    -- create a new Regent variable named resulti (where i is the loop variable, 
    -- e.g. ‘result0’) and store it in table ‘results’
    results[i] = regentlib.newsymbol(double, "result"..i)
    -- set result values to zero in preparation for accumulation
    statements:insert(rquote var [results[i]] = 0.0 end)
  end

  local pattern12 = generateJFockSpinPattern(L12)
  local pattern34 = generateJFockSpinPattern(L34)
  for q = 0, H34-1 do -- inclusive
    for p = 0, H12-1 do -- inclusive
      local Np, Lp, Mp = unpack(pattern12[p+1])
      local Nq, Lq, Mq = unpack(pattern34[q+1])
      local N, L, M = Np + Nq, Lp + Lq, Mp + Mq
      if (Nq + Lq + Mq) % 2 == 0 then
        statements:insert(rquote
          [results[p]] += density[q] * [R[N][L][M][0]]
        end)
      else
        statements:insert(rquote
          [results[p]] -= density[q] * [R[N][L][M][0]]
        end)
      end
    end
  end
  for i = 0, H12-1 do -- inclusive
    statements:insert(rquote
      accumulator[i] += [results[i]]
    end)
  end
  return statements
end
