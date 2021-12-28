import "regent"

local rsqrt = regentlib.rsqrt(double)
local round = regentlib.llrint(double)
local exp = regentlib.exp(double)
local SQRT_PI = math.sqrt(math.pi)

-----------------------------------------------------------------------------------------
-- Terra function that generates statements that iterpolate the Boys function to compute 
-- base auxiliary values R000J
-----------------------------------------------------------------------------------------
function generateStatementsComputeBoys(boys, length, t, r_gamma_table)

  local factors = {
    [0] = 2.7166345615120040E+01,
    [1] = 2.6273691011184581E+01,
    [2] = 2.5532392190003377E+01,
    [3] = 2.4901182676949695E+01,
    [4] = 2.4353368394557731E+01,
    [5] = 2.3870740944052123E+01,
    [6] = 2.3440355360724957E+01,
    [7] = 2.3052697787043194E+01,
    [8] = 2.2700582728555508E+01,
    [9] = 2.2378458056330899E+01,
    [10] = 2.2081949880679126E+01,
    [11] = 2.1807554706399655E+01,
    [12] = 2.1552425359950535E+01,
    [13] = 2.1314218509958917E+01,
    [14] = 2.1090983755524045E+01,
    [15] = 2.0881081442144000E+01,
    [16] = 2.0683120753173197E+01,
    [17] = 0.0000000000000000E+00,
  }
  local downwardsRecursion = terralib.newlist({rquote
    var t_scaled : double = [factors[length-1]] * t
    var coefficients : double[5] = r_gamma_table[{length-1, round(t_scaled)}];
    [boys[length-1]] = coefficients[0] + t_scaled * (coefficients[1]
                                       + t_scaled * (coefficients[2]
                                       + t_scaled * (coefficients[3]
                                       + t_scaled * coefficients[4])))
  end})
  for j = length-2, 0, -1 do -- inclusive
    downwardsRecursion:insert(rquote
      [boys[j]] = (2.0 * t * [boys[j+1]] + exp(-t)) * (1.0 / (2 * j + 1))
    end)
  end

  local upwardsRecursionAsymptotic = terralib.newlist({rquote
    [boys[0]] = rsqrt(t) * (SQRT_PI / 2.0)
  end})
  for j = 0, length-2 do -- inclusive
    upwardsRecursionAsymptotic:insert(rquote
      [boys[j+1]] = ((2 * j + 1) / 2.0) * [boys[j]] * (1.0 / t)
    end)
  end

  local statements = terralib.newlist()
  for j = 0, length-1 do -- inclusive
    statements:insert(rquote var [boys[j]] end)
  end
  statements:insert(rquote
    if t < 25 then
      [downwardsRecursion]
    else
      [upwardsRecursionAsymptotic]
    end
  end)
  return statements
end

---------------------------------------------------------------------------------------------
-- Terra function that generates statements that compute Hermite polynomials given by:
--    R000J = lambda * (-2*alpha)^J * F_J(t)
--    R00MJ = c * R00(M-1)(J+1) + (M-1) * R00(M-2)(J+1)
--    R0LMJ = b * R0(L-1)M(J+1) + (L-1) * R0(L-2)M(J+1)
--    RNLMJ = a * R(N-1)LM(J+1) + (N-1) * R(N-2)LM(J+1)
-- where F_J(t) is the Boys function. (See eqs. 14 in main text)
--
-- For example, an input length = 3 returns statements that interpolate the Boys fn. and 
-- compute R000J for J=1,2 (not shown) plus the following statement list:
--    R[0][0][1][0] = c * R[0][0][0][1]
--    R[0][0][1][1] = c * R[0][0][0][2]
--    R[0][0][2][0] = c * R[0][0][1][1] + 1.0 * R[0][0][0][1]
--    R[0][1][0][0] = b * R[0][0][0][1]
--    R[0][1][0][1] = b * R[0][0][0][2]
--    R[0][1][1][0] = b * R[0][0][1][1]
--    R[0][2][0][0] = b * R[0][1][0][1] + 1.0 * R[0][0][0][1]
--    R[1][0][0][0] = a * R[0][0][0][1]
--    R[1][0][0][1] = a * R[0][0][0][2]
--    R[1][0][1][0] = a * R[0][0][1][1]
--    R[1][1][0][0] = a * R[0][1][0][1]
--    R[2][0][0][0] = a * R[1][0][0][1] + 1.0 * R[0][0][0][1]
---------------------------------------------------------------------------------------------
function generateStatementsComputeRTable(R, length, t, alpha, lambda, a, b, c, r_gamma_table)
  -- interpolate Boys function and compute R000j values
  local statements = generateStatementsComputeBoys(R[0][0][0], length, t, r_gamma_table)
  for j = 0, length-1 do -- inclusive
    statements:insert(rquote
      [R[0][0][0][j]] *= lambda
      lambda *= -2.0 * alpha;
    end)
  end
  -- compute RNLM0 values from R000j base values 
  for N = 0, length-1 do -- inclusive
    for L = 0, length-1-N do -- inclusive
      for M = 0, length-1-N-L do -- inclusive
        for j = 0, length-1-N-L-M do -- inclusive
          if N == 0 and L == 0 and M == 0 then
            -- Do nothing. R000j has already been computed.
          elseif N == 0 and L == 0 and M == 1 then
            statements:insert(rquote
              var [R[0][0][1][j]] = c * [R[0][0][0][j+1]]
            end)
          elseif N == 0 and L == 0 then
            statements:insert(rquote
              var [R[0][0][M][j]] = c * [R[0][0][M-1][j+1]] + (M-1) * [R[0][0][M-2][j+1]]
            end)
          elseif N == 0 and L == 1 then
            statements:insert(rquote
              var [R[0][1][M][j]] = b * [R[0][0][M][j+1]]
            end)
          elseif N == 0 then
            statements:insert(rquote
              var [R[0][L][M][j]] = b * [R[0][L-1][M][j+1]] + (L-1) * [R[0][L-2][M][j+1]]
            end)
          elseif N == 1 then
            statements:insert(rquote
              var [R[1][L][M][j]] = a * [R[0][L][M][j+1]]
            end)
          else
            statements:insert(rquote
              var [R[N][L][M][j]] = a * [R[N-1][L][M][j+1]] + (N-1) * [R[N-2][L][M][j+1]]
            end)
          end
        end
      end
    end
  end

  return statements
end
