#!/usr/bin/env python

import numpy as np
import os
from astropy.io import fits

# Read FITS start times
def get_starttimes(froot):
    """Retrieve information about PSRFITS files

    Input:
       froot: Root of PSRFITS filenames, string
              e.g. fits/puppi_57614_C0531+33_0803 for fits/puppi_57614_C0531+33_0803_*.fits

    Output:
       files: array of filenames in PSRFITS observation
       stt_imjd: numpy array of integer MJDs of observation start
       stt_smjd: numpy array of integer seconds of observation start
       stt_offs: numpy array of offset from start time (seconds)
       nsblk: numpy array of bumber of spectra per subint
       tbin: numpy array of sampling time (seconds)
       nsuboffs: numpy array of subint offset with respect to observation start
       nsub: numpy array of number of subints within each file

    """  
    # Start file
    ifile=1

    # Empty lists
    files=[]
    stt_imjd=[]
    stt_smjd=[]
    stt_offs=[]
    tbin=[]
    nsblk=[]
    nsub=[]
    nsuboffs=[]

    # Loop over files
    while True:
        # Format filename
        fname="%s_%04d.fits"%(froot,ifile)

        # Check if the file exists
        if os.path.isfile(fname):
            # Increment counter
            ifile+=1

            # Store file name
            files.append(fname)

            # Read headers
            fits_file=fits.open(fname,memmap=True)
            fits_hdr=fits_file[0].header
            subint_hdr=fits_file['SUBINT'].header                                                                                              

            # Get information
            stt_imjd.append(fits_hdr['STT_IMJD'])
            stt_smjd.append(fits_hdr['STT_SMJD'])
            stt_offs.append(fits_hdr['STT_OFFS'])
            tbin.append(subint_hdr['TBIN'])
            nsblk.append(subint_hdr['NSBLK'])
            nsuboffs.append(subint_hdr['NSUBOFFS'])
            nsub.append(subint_hdr['NAXIS2'])
        else:
            break

    return files,np.array(stt_imjd),np.array(stt_smjd),np.array(stt_offs),np.array(nsblk),np.array(tbin),np.asarray(nsuboffs),np.asarray(nsub)

# Extract subints from a single file
def extract_subints_from_single_file(infname,outfname,isubmin,isubmax):
    """Extract subints from a single FITS file and store them as a different FITS file

    Input:
       infname: Name of input FITS file
       outfname: Name of output FITS file
       isubmin,isubmax: subint range to extract and store

    """  
    # Open file
    fits_file=fits.open(infname,memmap=True)

    # Read HDU 0 and 1 and their headers
    fits_hdu=fits_file[0]
    fits_hdr=fits_hdu.header
    subint_hdu=fits_file[1]
    subint_hdr=subint_hdu.header

    # Copy subint header
    new_subint_hdr=subint_hdr

    # Adjust NSUBOFFS
    new_subint_hdr['NSUBOFFS']+=isubmin

    # Create a new primary HDU and binary table HDU
    new_fits_hdu=fits.PrimaryHDU(data=None,header=fits_hdr)
    new_subint_hdu=fits.BinTableHDU(data=subint_hdu.data[isubmin:isubmax],header=new_subint_hdr)
    new_fits_file=fits.HDUList([new_fits_hdu,new_subint_hdu])

    # Write to new file
    new_fits_file.writeto(outfname,clobber=True)

    # Close files
    fits_file.close()
    new_fits_file.close()

    return


# Extract subints from an observation
def extract_subints_from_observation(froot,path,tbursts,isub0,isub1):
    """Extract subints from a PSRFITS observation

    Input:
       froot: root file PSRFITS filename
       path: output pat
       tbursts: list of burst times to extract
       isub0: number of subints to extract before the pulse (negative)
       isub1: number of subints to extract after the pulse

    """  
    # Get start times
    files,stt_imjd,stt_smjd,stt_offs,nsblk,tbin,nsuboffs,nsub=get_starttimes(froot)

    # File offsets
    offsets=nsblk*nsuboffs*tbin

    # Start and end times
    tstart=offsets
    tend=offsets+nsblk*nsub*tbin

    # Loop over bursts
    for tburst in tbursts:
        # Skip bursts outside of observation
        if tburst<0.0 or tburst>np.max(tend):
            print "Burst time %f not in range of observation!"%tburst
            continue

        # Loop over files
        for i in xrange(len(files)):
            # Time in between two files
            if tburst>=tstart[i] and tburst<tend[i]:
                toff=tburst-tstart[i]
                isub=int(np.floor(toff/(nsblk[i]*tbin[i])))
                
                # Subint limits
                isubmin=isub+isub0
                isubmax=isub+isub1

                # Some logic for dealing with file breaks
                if isubmin>=0 and isubmax<nsub[i]:
                    # Output filename
                    fname="%s_%08.3f.fits"%(os.path.join(path,os.path.basename(froot)),tburst)
                    
                    print "Extracting subints %03d to %03d from %s to %s"%(isubmin,isubmax,files[i],fname) 
                    
                    # Extract subints
                    extract_subints_from_single_file(files[i],fname,isubmin,isubmax)
                else:
                    print "Pulse at t=%g,isub=%d extends over a file break."%(tburst,isub)
                        
    return
    
if __name__ == '__main__':     

    # File root name    
    froot="../../../puppi_57614_C0531+33_0803"

    path="extracted_pulses"

    # Subint offsets
    isub0=-2
    isub1=+8
    
    # Burst times
    #    tbursts=np.array([1970.861179,318.310236,2367.424102,6033.051320])
    tbursts=np.array([214.0])

    # Extract subints
    extract_subints_from_observation(froot,path,tbursts,isub0,isub1)
