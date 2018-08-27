import visa
import time
import datetime
import numpy as np
import csv
import pandas as pd
import time
import glob
import sys
import matplotlib.pyplot as plt


class TransientDetector():
    def __init__(self):
        self.rm = visa.ResourceManager()
        self.gpibConn = self.rm.open_resource(
            'GPIB0::23::INSTR', read_termination='\r\n')
        self.outputFileName = None
        self.outputFile = None
        self.headerStored = False
        self.IsFieldThresholdScanSetup = False

    def clearBuffer(self):
        # Read the rest of the buffer to clear things out
        buff_clear = 0
        read_data = ""
        while(buff_clear == 0):
            try:
                read_data = self.gpibConn.read_raw()
                print(("Clearing the buffer: %s" % read_data))
            except:
                #print("Buffer is clear. Last value %s"%(read_data))
                buff_clear = 1

    def setupLabThresholdScan(self, outputFileName, minThresholdmV=0.00, maxThresholdmV=300, deltaThresholdmV=1.0):
        print("setting up threshold scan")
        # store the settings of the threshold scan.
        self.minThresholdmV = minThresholdmV
        self.maxThresholdmV = maxThresholdmV
        self.deltaThresholdmV = deltaThresholdmV
        self.NThresholds = int(
            np.ceil(self.maxThresholdmV/self.deltaThresholdmV))
        self.countPeriod = 1  # hard coding this now. This is in seconds
        # number of periods to count, we're manually doing the scan, so set this to one
        self.NPeriods = 1

        # First read anything in the buffer to clear things out
        self.clearBuffer()

        # Set up the photon counter
        # Hard-coding the following settings
        # Assuming that the pulse to be counted is connected to input 1
        # Counting on counter A and the T counter counts the internal clock
        # The count window is set to 1s.
        # Assuming that the gate has been set up using the front panel if you want it.
        self.gpibConn.write("CM0")  # counting mode A,B for T Preset
        self.gpibConn.write("CI0,1")  # counter A to input 1
        self.gpibConn.write("CI2,0")  # counter T to 10 MHz
        self.gpibConn.write("CP2,1E7")  # set T Set to 1s
        # note: T set is the number of cycles of the 10MHz clock, not seconds
        self.gpibConn.write('NP'+str(self.NPeriods))  # set number of periods
        self.gpibConn.write('DM0,0')  # set discriminator mode to fixed

        # write the header information to the file
        self.outputFileName = outputFileName
        if(self.outputFileName == 'None'):
            ts = time.time()
            st = datetime.datetime.fromtimestamp(
                ts).strftime("%Y-%m-%d_%H-%M-%S")
            self.outputFileName = "photonCounter_%s.csv"

        print(("Writing to output file %s" % (self.outputFileName)))
        self.outputFile = open(self.outputFileName, 'w')
        self.outputFile.write("Threshold Scan (%1.3f, %1.3f, %1.3f) mV, " % (
            self.minThresholdmV, self.maxThresholdmV, self.deltaThresholdmV))
        self.outputFile.write("Count Period %1.3e s" % self.countPeriod)

    def runLabThresholdScan(self):

        if(self.outputFileName == None):
            self.setupLabThresholdScan()

        self.outputFile.write('Discriminator Level (V), Counts\n')
        print("Running scan")
        D = np.arange(self.minThresholdmV, self.maxThresholdmV,
                      self.deltaThresholdmV)
        self.gpibConn.write("DL0,0")  # set discriminator threshold to zero
        #read_data = self.gpibConn.read_raw()
        # print(read_data)

        for discValue in D:
            # Set the discriminator threshold and then read
            # the associated data. You'll need to read the
            # discriminator threshold back and wait for the
            # counter to be ready.

            # set the discriminator threshold
            # set the discriminator threshold to the first one
            self.gpibConn.write("DL0,"+str(discValue/1.e3))
            # every time you write to the discriminator level, have to read it back
            read_disc = self.gpibConn.query("DL0")
            # print(read_disc)

            # start the scan
            self.gpibConn.write("DCL")  # clear the device
            self.gpibConn.write('CR;CS')  # start the scan

            # Sleep for a few seconds. uncomment this and below to verify syncing
            # time.sleep(0.9)

            # Read the data status bit and wait for data to be ready
            ss_int = 0
            ss = self.gpibConn.query("SS1")
            while (ss_int != 1):
                try:
                    ss_int = int(ss)
                    ss = self.gpibConn.query("SS1")
                    ss_int = int(ss)
                    # print(ss_int)
                except ValueError:
                    print(("Error reading status bit: ", ss))
                    ss = self.gpibConn.read_raw()
                    ss_int = int(ss)
                    continue

            A = int(self.gpibConn.query("QA1"))
            if(A < 0):
                A = int(self.gpibConn.query("QA1"))

            print(('disc level: '+str(discValue)+' (mV),  counts: '+str(A)))
            self.outputFile.write(str(discValue)+','+str(A)+'\n')
            # Sleep for half a second. Uncomment to verify syncing
            # time.sleep(.5)

        self.outputFile.close()

    # Collect Header Information from User
    def getCommonInfo(self, run=None, freqband=None, pol=None,
                      weather=None, daynight=None, boxno=None, before_stage1_attenfilter=None, after_stage2_attenfilter=None, comments=None):
        # Store the settings of the run.
        # We want to store the start date and time
        # UHF/VHF boxes, antennas, polarization
        # Any comments on the run setup

        # time stamp
        self.timestamp = time.time()
        self.datestr = datetime.datetime.fromtimestamp(
            self.timestamp).strftime("%Y-%m-%d_%H-%M-%S")
        # store run number, if given
        # if not, then check this file folder and make sure you're not overwriting a run
        self.run = run
        # Add this later
        # if( self.run == None):
        #    runsInDirc = glob.glob(self.dirc + "run[0-9]+.*.hdf5")
        #    print runsInDirc
        if(self.run == None):
            self.run = int(eval(input("Run Number? (Integer please)")))

        # use the run number and date-time to name the output file
        self.outputFileName = "photonCounter_run%d_%s.hdf5" % (self.run, self.datestr)
        print(("File name %s" % self.outputFileName))

        # Now store the prompts
        # Frequency band -- either UHF or VHF
        self.freqband = freqband
        if((self.freqband == None) or (self.freqband != 'UHF' and self.freqband != 'VHF')):
            self.freqband = input("Frequency band? (UHF/VHF)")

        # Sets of Boxes -- either 1 or 2
        self.boxno = boxno
        if((self.boxno == None) or (self.boxno != 1 and self.boxno!= 2)):
            self.boxno = int(float(eval(input("Box No.? (1/2)"))))

        # Polarization -- either H or V
        self.pol = pol
        if((self.pol == None) or (self.pol != 'H' and self.pol != 'V')):
            self.pol = input("Polarization? (H/V)")

        # Attenuation or filters
        # Before stage1
        self.before_stage1_attenfilter = before_stage1_attenfilter
        if (self.before_stage1_attenfilter == None):
            self.before_stage1_attenfilter = input(
                "Attenuation or Filters before Stage 1? ( e.g. 0 dB, SLP-80)")

        # After Stage2
        self.after_stage2_attenfilter = after_stage2_attenfilter
        if (self.after_stage2_attenfilter == None):
            self.after_stage2_attenfilter = input(
                "Attenuation or Filters after Stage 2? ( e.g. 10 dB, 450 MHz Notch)")

        # Weather -- whatever you like
        self.weather = weather
        if((self.weather == None)):
            self.weather = input("Weather?")
        # Day / Night
        self.daynight = daynight
        if(self.daynight == None):
            self.daynight = input("Day or Night? (Day/Night)")
        # Comments
        self.comments = comments
        if(self.comments == None):
            self.comments = input("Any other comments?")

    # Setup run at a Fixed Threshold
    def setupFixedRun(self, dirc, runTimeMinutes=1.0, discThresholdmV=None, run=None, freqband=None, pol=None,
                      weather=None, daynight=None, before_stage1_attenfilter=None,
                      after_stage2_attenfilter=None, boxno=None, comments=None):

        # directory where you will store the run
        self.dirc = dirc
        if(dirc[-1] != "/"):  # make sure there's a slash at the end
            self.dirc = self.dirc + "/"

        # Ask the user for things:
        self.getCommonInfo(run=run, freqband=freqband, pol=pol,
                           weather=weather, daynight=daynight,
                           before_stage1_attenfilter=before_stage1_attenfilter, after_stage2_attenfilter=after_stage2_attenfilter, boxno=boxno, comments=comments)

        # Discriminator Threshold
        self.discThresholdmV = discThresholdmV
        if(self.discThresholdmV == None):
            self.discThresholdmV = float(
                eval(input("Discriminator Threshold (mV)?")))

        # Before we start asking the photon counter things
        # Read anything in the buffer to clear things out
        self.clearBuffer()

        # Now set up the run on the photon counter
        # and store its settings

        # runTime
        self.runTimeMinutes = runTimeMinutes
        self.countPeriodSeconds = 1.
        self.NPeriods = int(self.runTimeMinutes * 60. /
                            self.countPeriodSeconds)

        # Set up the photon counter
        # Hard-coding the following settings
        # Assuming that the signal to be counted is connected to input 1
        # Counting on counter A and the T counter counts the internal clock
        # The count window is set to 1s.
        # Assuming that the gate has been set up using the front panel if you want it.
        self.gpibConn.write("CM0")  # counting mode A,B for T Preset
        self.countMode = self.gpibConn.query("CM")

        self.gpibConn.write("CI0,1")  # counter A to input 1
        self.gpibConn.write("CI2,0")  # counter T to 10 MHz
        self.countAInput = self.gpibConn.query("CI0")
        self.countTInput = self.gpibConn.query("CI2")

        # set T Set to 1s, count period is hard-coded!!
        self.gpibConn.write("CP2,1E7")
        self.countPeriodSeconds = float(self.gpibConn.query(
            "CP2"))/10e6  # divide by clock on trigger T input
        # note: T set is the number of cycles of the 10MHz clock, not seconds

        # self.gpibConn.write('NP'+str(self.NPeriods)) #set number of periods
        # fix the number of periods to 1 on the counter, but do a manual scan
        self.gpibConn.write("NP1")
        #self.NPeriods = self.gpibConn.query("NP")

        self.gpibConn.write('DM0,0')  # set discriminator mode to fixed
        self.discAMode = self.gpibConn.query("DM0")

        # set the discriminator value
        self.gpibConn.write("DL0,%1.3f" % (self.discThresholdmV/1.e3))
        self.fixedDiscriminatorThresholdmV = float(
            self.gpibConn.query("DL0"))*1e3

        # Store all the header information in a dataframe
        self.header = pd.DataFrame({'run': run, 'timestamp': self.timestamp, 'filename': self.outputFileName,
                                    'freqband': self.freqband, 'pol': self.pol, 'weather': self.weather,
                                    'daynight': self.daynight, 'comments': self.comments,
                                    'boxno': self.boxno, 'before_stage1_attenfilter': self.before_stage1_attenfilter, 'after_stage2_attenfilter': self.after_stage2_attenfilter,
                                    'fixedDiscThresholdmV': self.fixedDiscriminatorThresholdmV, 'discAMode': self.discAMode,
                                    'runTimeMinutes': self.runTimeMinutes, 'countPeriodSeconds': self.countPeriodSeconds, 'NPeriods': self.NPeriods,
                                    'countMode': self.countMode, 'counterAInput': self.countAInput, 'counterTInput': self.countTInput,
                                    }, index=[0])

        print(("Writing data to %s%s" % (self.dirc, self.outputFileName)))
        self.header.to_hdf(self.dirc + self.outputFileName,
                           key='header', mode='a', )
        self.headerStored = True

    # Run at a Fixed Thresold
    def runFixedCounter(self, dirc, runTimeMinutes=1.0):

        # Make sure to store the header information
        if(self.headerStored == False):
            self.setupFixedRun(dirc=dirc, runTimeMinutes=runTimeMinutes)

        # First read anything in the buffer to clear things out
        self.clearBuffer()

        # Set up arrays for runs
        # Note that the timestamps are coming from the computer so they
        # are not synced to a GPS time or anything
        # Uncertainty is related to the latecy between sending and receiving a command
        self.timestamps = []
        self.counts = []
	
        print("Running for %2.1f minutes at discriminator threshold %2.2f mV"%(runTimeMinutes, self.fixedDiscriminatorThresholdmV))
        # start the scan and read the data back
        for i in range(self.NPeriods):
            # Store the time stamp
            self.timestamps.append(time.time())
            self.gpibConn.write('CR;CS')  # start the scan

            # Read the data status bit and wait for data to be ready
            ss_int = 0
            ss = self.gpibConn.query("SS1")
            while (ss_int != 1):
                try:
                    ss_int = int(ss)
                    ss = self.gpibConn.query("SS1")
                    ss_int = int(ss)
                    # print(ss_int)
                except ValueError:
                    print(("Error reading status bit: ", ss))
                    ss = self.gpibConn.read_raw()
                    ss_int = int(ss)
                    continue

            A = int(self.gpibConn.query("QA1"))
            if(A < 0):
                A = int(self.gpibConn.query("QA1"))
            self.counts.append(A)

        # make an array with the dsicriminator values
        self.discriminatorThresholdmV = np.zeros(
            len(self.counts)) + self.fixedDiscriminatorThresholdmV

        # Store all the run information in a dataframe
        self.runCounts = pd.DataFrame({'rates': np.array(
            self.counts), 'timestamp': self.timestamps, 'thresholdmV': self.discriminatorThresholdmV})

        print(("Writing data to %s%s" % (self.dirc, self.outputFileName)))
        self.runCounts.to_hdf(
            self.dirc + self.outputFileName, key='rates', mode='a')

        # Reset flag in case you want to run again
        self.headerStored = False

    # Setup run at a Threshold Scan in the Field
    def setupFieldThresholdScan(self, dirc, minThresholdmV=0.00, maxThresholdmV=300.0, deltaThresholdmV=1.0,
                                run=None, freqband=None, pol=None,weather=None, daynight=None,
                                before_stage1_attenfilter=None,after_stage2_attenfilter=None, boxno=None, comments=None):

        # directory where you will store the run
        self.dirc = dirc
        if(dirc[-1] != "/"):  # make sure there's a slash at the end
            self.dirc = self.dirc + "/"

        # Ask the user for things:
        self.getCommonInfo(run=run, freqband=freqband, pol=pol,weather=weather, daynight=daynight,
                                before_stage1_attenfilter=before_stage1_attenfilter,after_stage2_attenfilter=after_stage2_attenfilter, boxno=boxno, comments=comments)

        # Before we start asking the photon counter things
        # Read anything in the buffer to clear things out
        self.clearBuffer()

        # Now set up the run on the photon counter
        # and store its settings

        # runTime
        self.countPeriodSeconds = 1.

        # store the settings of the threshold scan.
        self.minThresholdmV = minThresholdmV
        self.maxThresholdmV = maxThresholdmV
        self.deltaThresholdmV = deltaThresholdmV
        self.NThresholds = int(
            np.ceil(self.maxThresholdmV/self.deltaThresholdmV))
        self.countPeriod = 1  # hard coding this now. This is in seconds
        # number of periods to count, we're manually doing the scan, so set this to one
        self.NPeriods = 1

        # First read anything in the buffer to clear things out
        self.clearBuffer()

        # Set up the photon counter
        # Hard-coding the following settings
        # Assuming that the pulse to be counted is connected to input 1
        # Counting on counter A and the T counter counts the internal clock
        # The count window is set to 1s.
        # Assuming that the gate has been set up using the front panel if you want it.
        self.gpibConn.write("CM0")  # counting mode A,B for T Preset
        self.countMode = self.gpibConn.query("CM")
        self.gpibConn.write("CI0,1")  # counter A to input 1
        self.gpibConn.write("CI2,0")  # counter T to 10 MHz
        self.countAInput = self.gpibConn.query("CI0")
        self.countTInput = self.gpibConn.query("CI2")
        self.gpibConn.write("CP2,1E7")  # set T Set to 1s
        # note: T set is the number of cycles of the 10MHz clock, not seconds
        self.countPeriodSeconds = float(self.gpibConn.query(
            "CP2"))/10e6  # divide by clock on trigger T input
        self.gpibConn.write('NP'+str(self.NPeriods))  # set number of periods
        self.NPeriods = self.gpibConn.query("NP")
        self.gpibConn.write('DM0,0')  # set discriminator mode to fixed
        self.discAMode = self.gpibConn.query("DM0")

        # Store all the header information in a dataframe
        self.header = pd.DataFrame({'run': run, 'timestamp': self.timestamp, 'filename': self.outputFileName,
                                    'freqband': self.freqband, 'pol': self.pol, 'weather': self.weather,
                                    'daynight': self.daynight, 'comments': self.comments,
                                    'discAMode': self.discAMode, 'fixedDiscThresholdmV': -1,
                                    'countPeriodSeconds': self.countPeriodSeconds, 'NPeriods': self.NPeriods,
                                    'countMode': self.countMode, 'counterAInput': self.countAInput, 'counterTInput': self.countTInput,
                                    }, index=[0])

        print(("Writing header to %s%s" % (self.dirc, self.outputFileName)))
        self.header.to_hdf(self.dirc + self.outputFileName,
                           key='header', mode='a', )
        self.headerStored = True
        self.IsFieldThresholdScanSetup = True

    # Run a Threshold Scan in the Field
    def runFieldThresholdScan(self, dirc):

        if(self.IsFieldThresholdScanSetup == False):
            self.setupFieldThresholdScan(dirc)

        # Define some variables to store
        Dset = np.arange(self.minThresholdmV,
                         self.maxThresholdmV, self.deltaThresholdmV)

        self.discriminatorThresholdmV = np.zeros(len(Dset))
        self.rates = np.zeros(len(Dset))
        self.timestamps = []
        print((len(self.discriminatorThresholdmV)))

        self.gpibConn.write("DL0,0")  # set discriminator threshold to zero
        #read_data = self.gpibConn.read_raw()
        # print(read_data)

        for i, discValue in enumerate(Dset):
            # Set the discriminator threshold and then read
            # the associated data. You'll need to read the
            # discriminator threshold back and wait for the
            # counter to be ready.

            # set the discriminator threshold
            # set the discriminator threshold to the first one
            self.gpibConn.write("DL0,"+str(discValue/1.e3))
            # every time you write to the discriminator level, have to read it back
            self.discriminatorThresholdmV[i] = self.gpibConn.query("DL0")

            # Store the time stamp
            self.timestamps.append(time.time())
            self.gpibConn.write('CR;CS')  # start the scan

            # Read the data status bit and wait for data to be ready
            ss_int = 0
            ss = self.gpibConn.query("SS1")
            while (ss_int != 1):
                try:
                    ss_int = int(ss)
                    ss = self.gpibConn.query("SS1")
                    ss_int = int(ss)
                    # print(ss_int)
                except ValueError:
                    print(("Error reading status bit: ", ss))
                    ss = self.gpibConn.read_raw()
                    ss_int = int(ss)
                    continue

            A = int(self.gpibConn.query("QA1"))
            if(A < 0):
                A = int(self.gpibConn.query("QA1"))
            self.rates[i] = A

        # Store all the run information in a dataframe
        self.runCounts = pd.DataFrame(
            {'timestamp': self.timestamps, 'rates': self.rates, 'thresholdmV': self.discriminatorThresholdmV})

        print(("Writing data to %s%s" % (self.dirc, self.outputFileName)))
        self.runCounts.to_hdf(
            self.dirc + self.outputFileName, key='rates', mode='a')

        # Reset flag in case you want to run again
        self.headerStored = False
        self.IsFieldThresholdScanSetup = False


if __name__ == "__main__":

    ########################################
    # Run settings

    # Run Mode
    run_mode = ''
    if len(sys.argv) > 1:
        run_mode = sys.argv[1]

    # Run Minutes
    run_minutes = 10.0
    if len(sys.argv) > 2:
        run_minutes = float(sys.argv[2])

    # path for the outputfile
    dircPrefix = "/home/radio/data/transDet/transientDetector/photonCounter"

    ########################################
    sr = TransientDetector()
    #sr.setupLabThresholdScan("/home/radio/data/transientDetector/photonCounter/May27_VHFBoxes1_VPOL_ThresholdScan_Terminated.csv", minThresholdmV=0.0, maxThresholdmV=100.0)
    # sr.runLabThresholdScan()
    if(run_mode == 'fixed_threshold'):
        sr.setupFixedRun(dirc=dircPrefix)
        sr.runFixedCounter(dirc=dircPrefix, runTimeMinutes=run_minutes)
    elif(run_mode == 'threshold_scan'):
        sr.runFieldThresholdScan(dirc=dircPrefix)
    else:
        print(
            "\tusage:ipython transientDetector.py threshold_scan\n \t\t--or--\n \t ipython transientDetector.py fixed_threshold runMinutes[default 10.0]\n")
