function setup_nctoolbox_java
%addjars(fileparts(which(mfilename)));
addjars('C:\Program Files\MATLAB\R2023a\java\jar');

% configure log4j. Log4j is included in Matlab's classpath
root = org.apache.log4j.Logger.getRootLogger();

% Create a ConsoleAppender object with the specified layout
layout = org.apache.log4j.PatternLayout('%d{ISO8601} [%t] %-5p %c %x - %m%n');
consoleAppender = org.apache.log4j.ConsoleAppender();
consoleAppender.setLayout(layout);

% Add the ConsoleAppender to the root logger
root.addAppender(consoleAppender);

%root.addAppender(org.apache.log4j.ConsoleAppender(org.apache.log4j.PatternLayout('%d{ISO8601} [%t] %-5p %c %x - %m%n')));
root.setLevel(org.apache.log4j.Level.WARN)