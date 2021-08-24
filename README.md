# DENA-overflights

The DENA-overflights repository includes two components: overflight audiblity ground-truthing and creation and analysis of a representative 'active space' corresponding to propeller noise for a geographic region. 

DENA-overflights has several dependencies, the only of which included here is `gdal_transform.exe` and its supporting files. Remaining dependencies are found at the Github links below: 
  > Iyore: <a href="https://github.com/nationalparkservice/iyore">NPS/Iyore</a> <br>
  > NMSIM-Python: <a href="https://github.com/dbetchkal/NMSIM-Python">dbetchkal/NMSIM-Python</a> <br>
  > soundDB: <a href="https://github.com/gjoseph92/soundDB">gjoseph92/soundDB</a>

Additionally, overflights are queried from a databased developed by <a href="https://github.com/smHooper/flightsdb">smHooper</a>. 

The DENA_overflights repository holds two function libraries `active_space_utils.py` and `noise_metrics.py`. 


## Ground Truthing

Ground truthing involves loading existing flight tracks (either GPS or ADS-B) on top of an existing acoustic record from soundDB. The ground truthing user interface loads flights from a study area around microphone coordinates and plots the (sparse) tracks on a map with a highlighted region of (dense, spline-fit) points. On the right, a spectrogram from the acoustic record shows the noise for the same temporal bounds. Moving the sliders allows the user to annotate when a flight is audible (in the time domain) and _see_ the result in the spatial domain. Accounts for the speed of sound (assumed to be 343 m/s). 
![truthing](https://user-images.githubusercontent.com/8905274/130698861-3172948d-217b-4d4d-ad5e-7bd8b4a8c00c.png)

The ground truthing results in a saved annotations file that can be loaded in later to distinguish between audible and inaudible points. The goal is to create a shape in NMSim that best fits the 'audible' region without including too much of the 'inaudible' region; this region is dubbed the _active space_. 

Note that the active space is truly a three-dimensional space and not a two-dimensional plane; that is, the area of a planar slice of the active space is dependent on its elevation (or flight altitude). Additionally, the active space corresponds to a particular noise source or aircraft; the active space for a helicopter is not the same as the active space for a fixed-wing aircraft or a commercial jet. 

A sample inventory of ground truthing sites in Denali performed in 2021: 
![inventory](https://user-images.githubusercontent.com/8905274/130700132-c95aba5c-00e6-4707-8be3-fea9ba28ba90.png)


## Active Space Creation

There are three key components to active space creation: acoustic ambience levels, noise source file, and digital elevation model (DEM). NMSim is a sound propagation model that simulates the propagation of sound from a noise source over terrain; comparing NMSim results to the acoustic ambience (either third-octave bands _or_ broadband) will determine whether or not the flight is audible in a given location. This usage of NMSim is certainly <a href="https://nwtteis.com/portals/nwtteis/files/references/Ikelheimer_2004_NMSim_User_Manual.pdf">outside of its realm of intended use</a> but does an astoundingly good job regardless. The process for creating an active space looks a bit like the following: 
![algorithm](https://user-images.githubusercontent.com/8905274/130701270-c6c26f0a-ff7d-4812-95d3-e93640b00fe5.png)

Because broadband ambience is a viable comparison to detect audibility (although from 2021 ground truthing, it was determined that a noise source of around +23 dB should be used to compensate), park-wide or state-wide broadband models can be used to compute the active space for _any_ location given a DEM. 


## Noise Impact Analysis and Metrics

This final step is under development. The end goal is to quantify the noise impact of _real overflights_ given the theoretical active space with one grand assumption: if a flight is within the active space, then it is audible. For active spaces with landings (e.g. into basecamp on Denali or in the Ruth amphitheater), we will see that this is not a good assumption, but this is future work. 

Given a flight track, interpolating/densifying points resolves when a flight crosses into and out of an active space. The time inside an active space is therefore a noise interval, and a single track could theoretically pass in and out of an active space multiple times in one flight. The noise impacts can be derived most simply by counting the number of non-overlapping noise events (intervals), and more complicated metrics can be derived from the noise-free intervals (times where there are no noise intervals) or even factoring in the distance of the overflight (energy based metrics). 

Note that this method of quantifying noise impacts only uses flight tracks and computational effort... no microphones in the field! Here is an example from the Triple Lakes trail on a particularly noisy afternoon in June, 2021: 
![metrics](https://user-images.githubusercontent.com/8905274/130703201-efd366cc-5a04-4443-9a02-6a8ba702f490.png)
