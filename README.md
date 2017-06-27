# ThicknessEstimationProgram
Estimates the thickness of 2D materials using optical microscope images taken with narrow bandpass filters
Built in Python 2.7

Run EstimationProgram.py using Python 2.7

Spyder 2.7 from Anaconda can be downloaded from here:
https://www.continuum.io/downloads (make sure you download Python version 2.7!)
Run Spyder once that is downloaded
Then, open 'EstimationProgram.py'

If it is your first time opening the program you will need to press F6 or click 'Run Settings' at the top.
Select "Execute in a new dedicated Python Console" at the top
Select "Interact with the Python console after execution" at the bottom
Click run.

After you have done this once, you can just click 'Run' or press F5

Quick questions:
1) Why can't I choose a PDMS thickness?

The PDMS will be so thick (>0.1mm) that it will be out of the depth of field of the microscope. Light reflected from the back of the PDMS will not reach the microscope aperture, so the PDMS can be taken as semi-infinite.
In the program, if the user specifies 'PDMS' then the refractive index values of Si and SiO2 are just taken to be that of PDMS, with an arbitrary value of 300nm assigned to the SiO2 layer because it is irrelevant.
Hence, no graphs are available of any variable against PDMS thickness - nothing would change with PDMS thickness as long as the PDMS is not of the order of ~100nm, which to my knowledge, no PDMS in the lab is.

2) What if I want to use a different substrate from Silicon or PDMS?

Write in the refractive index values as you would with another material. Then either tweak the 'Silicon or PDMS?' question for each option or just change the refractive index of PDMS ('nvalsPDMS') to the refractive index of your material. This is assuming you're using something like quartz or glass which will be on the order of >0.1mm, i.e. much thicker than the SiO2 layer on Si in the silicon substrate.

3) How do I add new materials?

This will be detailed shortly. While it is fairly trivial, unfortunately it is laborious.

4) X or Y is not working...

Firstly, try messing around with the program. Is it breaking down when it's asking you to enter a thickness? If so, have you tried to enter something such as '300' instead of '300nm'? Etc.
If it is not working and you are truly stuck, submit a pull request or comment to this and I may see it and get back to you.
I have tested the program thoroughly in its current form and everything should be working.

5) How can I get the contrast for heterostructures?

Currently you cannot, until I add a piece of the program that does the intensity and contrast calculation in the form of matrices. Hopefully I will get round to this soon.
