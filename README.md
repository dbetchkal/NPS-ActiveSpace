[![DOI](https://zenodo.org/badge/389775527.svg)](https://zenodo.org/badge/latestdoi/389775527)
# NPS-ActiveSpace

An ***active space*** is a well-known sensory concept from bioacoustics ([Marten and Marler 1977](https://www.jstor.org/stable/pdf/4599136.pdf), [Gabriele et al. 2018](https://www.frontiersin.org/articles/10.3389/fmars.2018.00270/full)). It represents a geographic volume whose radii correspond to the limit of audibility for a specific signal in each direction. In other words, an active space provides an answer to the question, *"how far can you hear a certain sound source from a specific location on the Earth's surface?"*

This repository is designed to estimate active spaces for motorized noise sources transiting the U.S. National Park System. Aircraft are powerful noise sources audible over vast areas. Thus [considerable NPS management efforts have focused on protecting natural quietude from aviation noise intrusions](https://www.nps.gov/subjects/sound/overflights.htm). For coastal parks, vessels are similarly powerful noise sources of concern. For both transportation modalities `NPS-ActiveSpace` provides meaningful, quantitative spatial guides for noise mitigation and subsequent monitoring. 

## Example

Consider an example active space, below. It was computed using data from a long term acoustic monitoring site in Denali National Park, DENAUWBT Upper West Branch Toklat ([Withers 2012](https://irma.nps.gov/DataStore/Reference/Profile/2184396)), indicated by the black point in the center of the image. The bold black polygon delineates an active space estimate for flights at 3000 meters altitude. Points interior to the polygon are predicted to be audible, those exterior, inaudible. <br> 

Superposed over the polygon are flight track points colored by their corresponding audibility status. `NPS-ActiveSpace` includes an application that leverages the acoustic record to ground-truth audibility of co-variate vehicle tracks from GPS databases. Ground-truthing is used to "tune" an active space to the appropriate geographic extent via mathematical optimization.<br>
<br>
<img src="https://github.com/dbetchkal/NPS-ActiveSpace/blob/main/nps_active_space/img/NPS-ActiveSpace_example.png" alt="active space polygon example" width="300">


## Packages

This project is made up of four packages:

[`utils`](https://github.com/dbetchkal/NPS-ActiveSpace/tree/main/nps_active_space#utils): diverse utilities - file I/O, geoprocessing computations, acoustic propagation modelling, and detection statistics
    
[`ground-truthing`](https://github.com/dbetchkal/NPS-ActiveSpace/blob/main/_DENA/README.md#ground-truthing): a `tkinter`-based ground-truthing application

[`active-space`](https://github.com/dbetchkal/NPS-ActiveSpace/blob/main/_DENA/README.md#generate-active-space): generate and optimize active space polygons

`analysis`: estimate acoustic metrics from the intersection of an active space polygon and vehicle tracks (*currently in development, see `Analysis` branch*)

Also included are noise source [data](https://github.com/dbetchkal/NPS-ActiveSpace/tree/v2/nps_active_space/data) for tuning active space polygons.

***For more specific information on each package, view their individual READMEs.***

## Order of Operations

While each package can be used and run individually, the project was designed so that outputs of one package work seamlessly as the inputs for another. 

Packages were designed to be run in the following order:

`ground-truthing` $\rightarrow$ `active-space` $\rightarrow$ `analysis`

---

## License

### Public domain

This project is in the worldwide [public domain](LICENSE.md):

> This project is in the public domain within the United States,
> and copyright and related rights in the work worldwide are waived through the
> [CC0 1.0 Universal public domain dedication](https://creativecommons.org/publicdomain/zero/1.0/).
>
> All contributions to this project will be released under the CC0 dedication.
> By submitting a pull request, you are agreeing to comply with this waiver of copyright interest.

## Publications

Publications about `NPS-ActiveSpace`:

>Betchkal, D.H., J.A. Beeco, S.J. Anderson, B.A. Peterson, and D. Joyce. 2023. Using Aircraft Tracking Data to Estimate the Geographic Scope of Noise Impacts from Low-Level Overflights Above Parks and Protected Areas. Journal of Environmental Management 348(15): 119201 https://doi.org/10.1016/j.jenvman.2023.119201
