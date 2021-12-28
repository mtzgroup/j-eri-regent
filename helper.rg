-----------------------------------------------------------------
-- Returns the maximum momentum determined at Regent compile time
-----------------------------------------------------------------
local _max_momentum = nil
function getCompiledMaxMomentum()
  if _max_momentum == nil then
    for i, arg_value in ipairs(arg) do
      if arg[i] == "-L" then
        _max_momentum = StrToL[arg[i+1]]
      end
    end
    assert(_max_momentum ~= nil,
           "Must give max angular momentum `-L [S|P|D|F|G]`")
  end
  return _max_momentum
end

-- The string representation of a given angular momentum
LToStr = {
  [0] = "S", [1] = "P", [2] = "D", [3] = "F", [4] = "G"
}
StrToL = {
  ["S"] = 0, ["P"] = 1, ["D"] = 2, ["F"] = 3, ["G"] = 4
}
LPairToStr = {
  [0] = "SS", [1] = "SP", [2] = "PP",
  [3] = "PD", [4] = "DD", [5] = "DF",
  [6] = "FF", [7] = "FG", [8] = "GG"
}

----------------------------------------------------------------
-- Returns the number of Hermite functions resulting from the
-- product of two Gaussians with angular momenta L1 + L2 = n
----------------------------------------------------------------
function tetrahedral_number(n)
  return n * (n + 1) * (n + 2) / 6
end

----------------------------------------------------------------
-- Returns a list of lists where `pattern[i] = {N, L, M}`.
--
-- For example, `generateJFockSpinPattern(2)` returns the table
-- 0 0 0
-- 1 0 0
-- 0 1 0
-- 0 0 1
-- 1 1 0
-- 1 0 1
-- 0 1 1
-- 2 0 0
-- 0 2 0
-- 0 0 2
-- Remember that indices of Lua lists start with 1.
----------------------------------------------------------------
function generateJFockSpinPattern(level)
  local pvec = {} -- create table
  for M = 0, level do -- inclusive
    for L = 0, level - M do -- inclusive
      for N = 0, level - M - L do -- inclusive
        table.insert(pvec, {N, L, M})
      end
    end
  end
  local pattern = {}
  for k = 0, level * level do -- inclusive
    for _, v in pairs(pvec) do -- iterator over pvec table
      if v[1] * v[1] + v[2] * v[2] + v[3] * v[3] == k then
        table.insert(pattern, v)
      end
    end
  end
  return pattern
end
