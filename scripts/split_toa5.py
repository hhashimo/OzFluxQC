'''
PURPOSE:
 Script to divide a Campbell Scientific TOA5 file into smaller chunks.
USE:
 At the command prompt in the subdirectory containing the scripty, type
  "python split_toa5_10hz.py"
 Edit the input file path, input file name, sampling period and output
 period to suit your needs.
 Typical values are sampling_period=0.1 (for 10 Hz) and output_period=1440
 for daily files.  EddyPro is able to read daily files but still produce
 output every 30 minutes.
USE CASE:
 The script was written to deal with TOA5 files of fast (10 Hz) data from the
 Wombat State Forest flux tower.  The raw TOB3 files from the CF card had
 been converted into TOA5 files (using the Loggernet CardConvert utility)
 without configuring CardConvert to write an output file every 30 minutes.
 The resulting TOA5 was 4.5 Gb in size and could not be input to EddyPro or
 easily handled by any other application (may have possible with the Loggernet
 utility Split).
AUTHOR: Peter Isaac (pisaac.ozflux@gmail.com)
DATE: 20th March 2014
'''
from datetime import datetime, timedelta

def get_datetime(line):
    # separate into parts using comma as delimiter
    parts = line.split(',')
    # strip the double quotes from the datetime string (first column)
    dt_str = parts[0].translate(None,'"')
    # convert the datetime string to a datetime object
    try:
        dt = datetime.strptime(dt_str,'%Y-%m-%d %H:%M:%S')
    except ValueError:
        dt = datetime.strptime(dt_str,'%Y-%m-%d %H:%M:%S.%f')
    # return the datetime object
    return dt

# things to change for each run
in_filepath = "../../Sites/Wombat/Data/Raw/2013/"
in_filename = "TOA5_3952.fast_std.dat"
sampling_period = 0.1         # sampling period, seconds eg 0.1 for 10 Hz, 0.05 for 20 Hz
output_period = 1440          # output period, minutes eg 30, 60 (for hourly), 1440 (for daily)

# start of main code
line = ''
sp_usecs = 0.1*1000000        # sampling period, microseconds
in_name = in_filepath + in_filename
infile = open(in_name,'rb')
# read the header lines
h1 = infile.readline()
h2 = infile.readline()
h3 = infile.readline()
h4 = infile.readline()
# read the first data line
line = infile.readline()
# get the datetime
dt = get_datetime(line)
# print the first datetime found in the file
print "First datetime in file is " + dt.strftime("%Y-%m-%d %H:%M:%S.%f")
# check to see if the file starts at the start of a day
if dt.hour==0 and dt.minute==0 and dt.second==0 and dt.microsecond==sp_usecs:
    # yay, file starts at start of day
    start_dt = dt                                           # get the start datetime of the current day
else:
    # boo, file doesn't start at start of day
    # get the start of the next day
    start_dt = datetime(dt.year,dt.month,dt.day,0,0,0)+timedelta(days=1)+timedelta(microseconds=sp_usecs)
    # loop over reading lines until we get to the start of the next day
    while dt<start_dt:
        line = infile.readline()
        dt = get_datetime(line)
# get the end datetime
end_dt = start_dt + timedelta(minutes=output_period)    # define the end datetime of this file
# print start and end datetimes to the screen
print "Start datetime of first file will be "+start_dt.strftime("%Y-%m-%d %H:%M:%S.%f")
print "End datetime of first file will be "+end_dt.strftime("%Y-%m-%d %H:%M:%S.%f")
# make the output file name
out_filename = in_filepath+'TOA5_'+start_dt.strftime("%Y%m%d%H%M")+'.dat'
# tell the user which day we are doing
print "Doing " + start_dt.strftime("%Y-%m-%d %H:%M")
# open the output file
outfile = open(out_filename,'w')
# write the header lines to the output file
outfile.write(h1)
outfile.write(h2)
outfile.write(h3)
outfile.write(h4)
# write out the first line if it was read while getting to the start of the first day
if line: outfile.write(line)
# big loop over the lines in the input file
while True:         # loop forever, we have to break out of loop on end of file
    # read the next line
    line = infile.readline()
    # line will be empty when we reach the end of file
    if not line: break
    # get the datetime
    dt = get_datetime(line)
    # check to see if we are at the end datetime for this period
    if dt==end_dt:
        # finished this period, close the output file
        outfile.close()
        # update the start datetime
        start_dt = end_dt
        # calculate the new end datetime
        end_dt = start_dt + timedelta(minutes=output_period)
        # generate the new output filename
        out_filename = in_filepath+'TOA5_'+start_dt.strftime("%Y%m%d%H%M")+'.dat'
        # tell the user which day we are doing
        print "Doing " + start_dt.strftime("%Y-%m-%d %H:%M")
        # open the file
        outfile = open(out_filename,'w')
        # write the header lines
        outfile.write(h1)
        outfile.write(h2)
        outfile.write(h3)
        outfile.write(h4)
        # write the first line
        outfile.write(line)
    else:
        # not at end of day, write line to output file
        outfile.write(line)
# close the files
infile.close()
outfile.close()
print "Done"