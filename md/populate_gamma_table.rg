import "regent"

local gamma_header = terralib.includec("md/gamma_table.h", {"-I", "./"})

local terra getGammaTable(j : int, index : int) : double[5]
  return gamma_header.gamma_table[j][index]
end

task populateGammaTable(r_gamma_table : region(ispace(int2d, {18, 700}), double[5]))
where
  writes(r_gamma_table)
do
  for index in r_gamma_table.ispace do
    r_gamma_table[index] = getGammaTable(index.x, index.y)
  end
end
