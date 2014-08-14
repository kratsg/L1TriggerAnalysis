require 'optparse'
options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: config.rb [options]"

  options[:file] = 'config'
  opts.on('-o', '--output FILE', 'Output File') { |v| options[:file] = v }

  options[:datasets] = 'datasets'
  opts.on('-F','--input DATASETS', 'Datasets to use') { |v| options[:datasets] = v}

  opts.on('-h', '--help', 'Display this screen') do
    puts opts
    exit
  end
end.parse!

localPython = File.expand_path "../.local"
mainPython  = File.expand_path "main.py"

#check that paths exist
unless File.exist?(localPython) then
  puts "Cannot find `.local/`. It should be in the parent directory `..`"
  exit
end

unless File.exist?(mainPython) then
  puts "Cannot find `main.py`. It should be in the same directory as this file."
  exit
end

configStart = <<-eos
##{Time.now}
# This file was dynamically generated from make_config.rb
Universe = vanilla 
Executable = main.sh
#Requirements = (OpSysAndVer =?= "SL6") && ( IS_RCC_uchicago )
Requirements = (OpSysAndVer =?= "SL6") && ( IS_RCC_fresnostate =!= True )
Error = out/error.$(Cluster)-$(Process)
Output = out/output.$(Cluster)-$(Process)
Log = out/log.$(Cluster)-$(Process)
+ProjectName="atlas-org-uchicago"
should_transfer_files = YES
when_to_transfer_output = ON_Exit
transfer_output         = True
transfer_input_files    = #{localPython}, #{mainPython}
transfer_output_files   = data
eos

start = 0
stop = 10000
step = 2000

numJobs = 0

open(options[:file], 'w+') do |f|
  f.puts(configStart)
  IO.readlines(options[:datasets]).map{|l| l.chomp}.each do |file|
    (start...stop).step(step) do |n|
      f.puts("Arguments = $(Process) #{file} #{n} #{step}")
      f.puts("Queue 1")
      numJobs += 1
    end
  end
end

puts "Datasets:\t#{File.expand_path options[:datasets]}"
puts "Output File:\t#{File.expand_path options[:file]}"
puts "#{numJobs} jobs created"
