import "regent"

require "helper"

local c = regentlib.c
local assert = regentlib.assert

struct Config {
  parallelism          : int;
  input_directory      : int8[512];
  parameters_filename  : int8[512];
  output_filename      : int8[512];
  verify_filename      : int8[512];
  num_trials           : int;
}

terra Config:dump()
  c.printf("\nConfig:")
  c.printf("\tparallelism: %d\n", self.parallelism)
  c.printf("\tinput_directory: %s\n", self.input_directory)
  c.printf("\toutput_filename: %s\n", self.output_filename)
  c.printf("\tverify_filename: %s\n", self.verify_filename)
  c.printf("\tnum_trials: %d\n\n", self.num_trials)
end

terra print_usage_and_abort()
  c.printf("Usage: regent coulomb.rg [OPTIONS]\n")
  c.printf("OPTIONS\n")
  c.printf("  -h                   : Print the usage and exit.\n")
  c.printf("  -L {value}           : Use {value} as the max momentum (set at compile time).\n")
  c.printf("  -i {directory}       : Use input files in {directory}.\n")
  c.printf("  -p {value}           : Partition bras {value} ways.\n")
  c.printf("  -o {file}            : Write output to {file}\n")
  c.printf("  -v {file}            : Verify the output with {file}\n")
  c.printf("  --trials {value}     : Run {value} times.\n")
  c.exit(0)
end

terra Config:initialize_from_command()
  self.parallelism = 1
  self.input_directory[0] = 0
  self.output_filename[0] = 0
  self.verify_filename[0] = 0
  self.num_trials = 1

  var args = c.legion_runtime_get_input_args()
  var i = 1
  while i < args.argc do
    if c.strcmp(args.argv[i], "-h") == 0 then
      print_usage_and_abort()
    elseif c.strcmp(args.argv[i], "-i") == 0 then
      if self.input_directory[0] ~= 0 then
        c.printf("Error: Only accepts one input directory!\n")
        c.abort()
      end
      i = i + 1
      c.strncpy(self.input_directory, args.argv[i], 512)
    elseif c.strcmp(args.argv[i], "-o") == 0 then
      if self.output_filename[0] ~= 0 then
        c.printf("Error: Only accepts one output file!\n")
        c.abort()
      end
      i = i + 1
      c.strncpy(self.output_filename, args.argv[i], 512)
    elseif c.strcmp(args.argv[i], "-v") == 0 then
      if self.verify_filename[0] ~= 0 then
        c.printf("Error: Only accepts one verify file!\n")
        c.abort()
      end
      i = i + 1
      c.strncpy(self.verify_filename, args.argv[i], 512)
    elseif c.strcmp(args.argv[i], "-p") == 0 then
      i = i + 1
      self.parallelism = c.atoi(args.argv[i])
    elseif c.strcmp(args.argv[i], "--trials") == 0 then
      i = i + 1
      self.num_trials = c.atoi(args.argv[i])
      assert(self.num_trials > 0, "Number of trials must be positive.")
    elseif c.strcmp(args.argv[i], "-L") == 0 then
      i = i + 1
    elseif c.strncmp(args.argv[i], "-ll:", 4) == 0 then
      c.printf("Passing flag %s to Legion runtime\n", args.argv[i])
      i = i + 1
    else
      c.printf("Warning: Unknown option %s\n", args.argv[i])
    end
    i = i + 1
  end

  if self.input_directory[0] == 0 then
    c.printf("Error: Input file not given!\n")
    print_usage_and_abort()
  end

  c.sprintf(self.parameters_filename, "%s/parameters.dat", self.input_directory)
end

return Config
