require 'logger'

module Logging

  @out = STDOUT
  @classname = ""
  @loggers = {}

  class << self

    def configure(config) # should be a log path, like /tmp/log.txt
      @out = config[:out] if (config[:out] and config[:out] != 'STDOUT')
      @classname = config[:class] if (config[:class])
    end

  end

  # This is the magical bit that gets mixed into your classes
  def log
    Logging.log
  end


  # Global, memoized, lazy initialized instance of a logger
  def self.log
    @log ||= Logger.new(@out)
    @log.level = Logger::INFO
    @log
    #@log.progname = @classname
    #@log
  end

end