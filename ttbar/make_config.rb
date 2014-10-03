require 'optparse'
require 'yaml'

options = {}
OptionParser.new do |opts|
  opts.banner = "Usage: config.rb [options]"

  options[:outputFile] = 'config'
  opts.on('-o', '--output FILE', 'Output File') { |v| options[:outputFile] = v }

  options[:configFile] = 'config.yml'
  opts.on('-F','--input CONFIG', 'Configuration options for me') { |v| options[:configFile] = v}

  opts.on('-h', '--help', 'Display this screen') do
    puts opts
    exit
  end
end.parse!

localPython = File.expand_path "../local.tar.gz"
mainPython  = File.expand_path "main.py"
x509Proxy = File.expand_path "/tmp/x509up_u#{%x(id -u).chomp}"

#check that paths exist
unless File.exist?(localPython) then
  puts "Cannot find `.local/`. It should be in the parent directory `..`"
  exit
end

puts "Found localPython\n\t#{localPython}"

unless File.exist?(mainPython) then
  puts "Cannot find `main.py`. It should be in the same directory as this file."
  exit
end

puts "Found main python script\n\t#{mainPython}"

unless File.exist?(x509Proxy) then
  puts "Cannot find your proxy file #{x509Proxy}.\n\tAre you sure you set it up? Run `localSetupFAX && voms-proxy-init -voms atlas` to set it up."
  exit
end

puts "Found x509 proxy\n\t#{x509Proxy}"

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
transfer_input_files    = #{localPython}, #{mainPython}, #{x509Proxy}
transfer_output_files   = data
environment             = "X509_USER_PROXY_FILENAME=x509up_u33155"
eos

configurations = YAML.load_file(options[:configFile])

totalJobs = 0

open(options[:outputFile], 'w+') do |f|
  f.puts(configStart)
  configurations['datasets'].each do |dataset|
    file = dataset['file']
    numEvents = dataset.has_key?('numEvents') ? dataset['numEvents'] : configurations['numEvents']
    numJobs = dataset.has_key?('numJobs') ? dataset['numJobs'] : configurations['numJobs']
    eventsPerJob = (numEvents/numJobs.to_f).ceil # round up in case, handle the edge case
    (0...numJobs).each do |i|
      f.puts("Arguments = $(Process) #{configurations['prefix']}#{file} #{i*eventsPerJob} #{[(i+1)*eventsPerJob, numEvents].min - i*eventsPerJob}")
      f.puts("Queue 1")
      totalJobs += 1
    end
  end
end

puts "Config File:\t#{File.expand_path options[:configFile]}"
puts "Output File:\t#{File.expand_path options[:outputFile]}"
puts "#{totalJobs} jobs created"
