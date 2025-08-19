## Rationale and Objectives  
In this project we aim to characterize sensor-level oscillatory activity from lag-modulated, repeated presentation of high-caloric food, compared with positive and neutral stimuli that respectively possess and lack perceptual salience.  

## Assumptions and Hypotheses   
Given prior discoveries (Codispoti et al., 2023) it is assumed that salient stimuli will exhibit enhanced desynchronization across repeated presentations, indicating larger cortical excitability.  

## Data Description  
Each participant from the 42 participants in the study was shown image repetitions.   
There were 18 conditions in the experiment, and images were categorized into 3 main semantic categories (food, positive and neutral). Each image was shown twice, with varying lags between repetitions (short, medium and long). The participants' brain activity was recorded using an MEG 4D-Neuroimaging system.   
Recordings were sampled at a rate of 1017 Hz and online band-pass filtered from 1–400 Hz. Data was preprocessed as follows: trial definition, line noise and artifact removal, low-pass filtering at 90 Hz, and ICA. The data was epoched in the time range of -0.3-0.8 s relative to event/stimuli onset.   

The data used in the project is:  
* Epoched data in a mat file named “ICA_dataorig”.
* Evoked data per experimental condition in a mat file "don_ERF".  
For each of the 42 study participants.  

The data for two of the 42 participants and can be found in the following link, for pipeline implementation purposes: [DSNS](https://livebiuac-my.sharepoint.com/:f:/g/personal/elizabeth_vaisman_live_biu_ac_il/EmTGDL0frsxFvlgsm3hs5woBhhSwPVyjT6H3Ak81gwvxgg?e=DkC58C)   
_Result of the implementation is summarised in a matlab figure "resulting_fig.fig" and can be downloaded from the same link._

## Implementation

Download the data and create the following directory structure:  
├───Main _Folder    
&emsp;│   dsns_script.m    
&emsp;├───TFR_DSNS      
&emsp;├───sub_003      
&emsp;&emsp;│   ICA_dataorig.mat      
&emsp;&emsp;├───ERF_within     
&emsp;&emsp;&emsp;|   don_ERF.mat    
&emsp;&emsp;├───TFR_within    
&emsp;└───sub_004      
&emsp;&emsp;│   ICA_dataorig.mat    
&emsp;&emsp;├───ERF_within   
&emsp;&emsp;&emsp;|   don_ERF.mat  
&emsp;&emsp;├───TFR_within    


### Changes in dsns_script.m for frequency analysis implementation 
 
 Change the following variables according your own paths and subject folder names.
 * Path variables: path_ft, path_files, datapath, datafile, savefile.  
 * Subject folder names:    

       % Define subjects and trigger values  
       subjects = {};  
       ranges = [3:16, 18:19, 21:22, 25:29, 31:34, 36:37, 39:41, 43:45, 47, 49:54];  
       for i = ranges  
           subjects{end+1} = sprintf('sub_%03d', i);  
       end
   
### Frequency Analysis Implementation Steps  
  
1.	Run the initialization, comment out the "% Atlas for reference part" till end of initialization (this part is used for cluster analysis).

2. Calculate induces data for each subject using calculate_induced, which subtract the ERF of each condition from the corresponding trials in ICA_dataorig.mat.
    
3.	Extract baseline and post-stimulus activity and preform a fourier transform on the data, using the "calculate_freq" function.
   
    _Notes_:   
    * _Baseline and post-stimulus frequency analysis are done seperately by changing the time_interval = [min max] variable and the  last input in calculate_freq to true for baseline and false otherwise._
    * Make sure to change: freq_resul for different desired frequency resolution and the lower frequency border, and fs indicating recording sampling frequency.
  
4. Concatenate the analysed post-stimulus and baseline data of all subjects using the "concatenate_freqs " function. "datafile" and "savefile" variables change accordingly.
   
    _Note: make sure load_file and save_file passed to the function match whther both regard baseline or post-stimulus. load_file is the same file that was saved
   in "calculate_freq" function._
  
6.	Load baseline and post-stimulus concatenated data for normalization.
   
7.	Normalize the induced post-stimulus data using "extract_2D_power" function, take the average power across all channels per frequency and then get the average power for a specific frequency band specified as an input to the function.
   
    _Note: result should be used in repeated measures ANOVA._   

8.	Plot the average power spectrum across all subjects for specific condition groups.
    
    _Note: "condGroups" variable contains the indices of the conditions listed in "trigVal". Changing the "condGroups" variable requires also changing  " customNames" according to the conditions' groups._  

trigVal: 8 – "oddball"  
Food: 10 – food-short repetition 1, 12 – food-medium 1, 14– food-long 1, 20 - food-short 2, 22 – food-medium 2,   
24 – food-long 2.   
Same for Positive: 110 112 114 120 122 124  
Same for Neutral: 210 212 214 220 222 224  

## Submittables   
A pdf report is available alongside the respective directories summarizing the main results from the whole pipeline implementation on the 42 subjects' data mentioned in the report.   
