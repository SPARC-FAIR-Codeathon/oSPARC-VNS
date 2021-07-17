
User usage:

## This should all be in the same wireframe (same /url): 

- select a device family using a drop-down menu
- select a particular device using a drop-down menu

- customise the selected device using input number fields : number of electrodes, [1]electrode size [3], electrode location [3], (which electrode is being edited, multiselect dropdown)
- set the size and shape of the insulator (electrode carrier)
- see a 2d graphical description of their selected device in 2 planes


- enter a URL to a nerve image (as MBF-XML, from sparc.science) to simulate
- upload a nerve image of their own (as MBF-XML) to simulate
- move / rotate the entered nerve cross-section in 2D so it sits inside the device
- zoom in and out in this view to see more or less context 

- choose what kind of simulation to run: {mesh only, fields only, generate axons only, simulate nerve recording}

- set the kind of axon population to use from a drop-down: {rat pelvic nerve, rat cervical vagus, human cervical vagus} 
(- eventually enter a URL or upload a file for axon population data)

- set what kind of nerve recording to simulate, in terms of axon activity levels {quasi-stationary, steadily changing, bursting}
(- set an input spike pattern to simulate by uploading a nwb or some other file?)

- get more information about their selected axon population (link)
- supply an email address to recieve emailed results

- Run the model: this will involve API calls to oSPARC plus the server keeping the session alive to know when to email or otherwise serve the data to the user 


## Visualisation might be under a seperate /url ? 

- visualise results (electrical fields, axon population, simulated nerve recording)stored in a file on the user's computer 
- set some basic properties about the data being visualised (drop-downs, multi-select drop-downs): which axon population is "on"? which electrode(s) are we looking at?  show waves, spikes, or both?
