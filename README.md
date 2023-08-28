# Success-Score

This repository provides a SOCCER©-independent way of creating Success-Scores based on position data. Please be aware that the code is far from package-ready and rather a simplistic proof of concept including major options / parameters of the Success-Score. With the implementation in Python the options for the Success-Score are now basically unlimited, allowing for a wide variety of visualizations and methodological adaptations. <br>

### The Success-Score
The Success-Score was developed back in 2017 and is one of many possible outputs provided by the match analysis software SOCCER© (Perl et al., 2013, Perl & Memmert, 2017). During my time at the Institute of Exercise Training and Sport Informatics at the German Sport University Cologne I had the opportunity to join the research project working with the Success-Score. The Success-Score essentially is intended to capture attacking (performance) in soccer by combining space and ball control using the concepts of efficiency and effort. For a detailed description of the Success-Score I recommend checking out the original paper (Perl & Memmert, 2017) and the first validation study (Brinkjans et al., 2022). Additional publications may follow. <br>

### The code
The code in this repository is divided into 3 Jupyer-notebooks:

+ *BallControl.ipynb* - Providing the code and methodology of the ball control component of the Success-Score, accompanied by some basic examples
+ *SpaceControl.ipynb* - Providing the code and methodology of the space control component of the Success-Score, accompanied by some basic examples
+ *SuccessScore.ipynb* - Providing the code and methodology of the Success-Score combining space and ball control, accompanied by some examples of the parameters available

In these scripts a work with a random data file of position data. Unfortunately I cannot provide the data. <br>

Additionally  a python file - *SuccessScore.py*, an adaptation of *floodlight's* space control code (Raabe et al., 2022, Source Code for Floodlight.Models.Space, 2023) is included. This file is essentially self-sufficient and includes all the necessary  functions to create Success-Scores. At this stage the code is dependent on position data and match information provided in a certain format (dfl - floodlight) and is build upon existing code from floodlight.Without but especially with the help of floodlight the adapation of the code to allow for the analysis of position data from different sources is easily implementable. Floodlight is a match analysis package developed at Institute of Exercise Training and Sport Informatics at the German Sport University Cologne which I highly recommend checking out: https://floodlight.readthedocs.io/en/latest/index.html | https://github.com/floodlight-sports/floodlight | Raabe et al., 2022 <br>

### References 
Brinkjans, D., Memmert, D., Imkamp, J., & Perl, J. (2022). Success-Score in Professional Soccer – Validation of a Dynamic Key Performance Indicator Combining Space Control and Ball Control within Goalscoring Opportunities. International Journal of Computer Science in Sport, 21(2), 32–42. https://doi.org/10.2478/ijcss-2022-0009 <br>
Perl, J., Grunz, A., Memmert, D., & Gutenberg-University, J. (2013). Tactics Analysis in Soccer – An Advanced Approach. International Journal of Computer Science in Sport, 12. <br>
Perl, J., & Memmert, D. (2017). A Pilot Study on Offensive Success in Soccer Based on Space and Ball Control – Key Performance Indicators and Key to Understand Game Dynamics. International Journal of Computer Science in Sport, 16(1), 65–75. https://doi.org/10.1515/ijcss-2017-0005 <br>
Raabe, D., Biermann, H., Bassek, M., Wohlan, M., Komitova, R., Rein, R., Groot, T. K., & Memmert, D. (2022). floodlight—A high-level, data-driven sports analytics framework. Journal of Open Source Software, 7(76), 4588. https://doi.org/10.21105/joss.04588 <br>
