## Project info

**Title:**
Integration of the single-trial time-resolved spectral connectivity (coherence; PLV) in Frites

**Project lead and collaborators:** 
Lima, Vinicius; Combrisson, Etienne (@kNearNeighbors)(@EtienneCmb)

**Description:**
The synchronization of the activity from distinct brain areas has been proposed to be one of the mechanisms by which them integrate while processing similar inputs in order to exchange information or encode the stimulus (Buzsaki G.,2006; Fries P., 2015). Based on this hypothesis, 
the time-course of the functional connectivity (dFC) can be measured from brain signals using metrics that capture their phase-relation such as the cross-spectra, 
the phase-locking value (PLV), and the coherence (Bastos A.M., Schoffelen J.M. 2006). Additionally, apart from estimating those metrics in a time-resolved 
manner, in order to be able to relate the dynamics of the phase-coupling and task-related behavioral events, it is also relevant to assess the dFC
at single-trial level, hence, avoiding averaging out non-phase-locked bursts of synchronization that are present in the dFC and may correspond to brain states relevant to determining, for instance, whether the information is being encoded during cognitive tasks by the coordinated activity of multiple cortical areas.  
Currently, xfrites - the testing repository associated to Frites (https://brainets.github.io/frites/) - has a function that estimates dFC in terms of the aforementioned metrics. 
For the present project, we aim to integrate it with FRITES and, more specifically, we aim to improve the documentation of the function, 
refine the current implementation and include code for unit testing. Other goals are to implement notebooks with examples that allow the user to have a better understanding of how to first, set the parameters to estimate the spectral connectivity and seconds, to interpret the metric's outcome, what are
its advantages and drawbacks.

**keywords:** Communication through coherence; spectral analysis; wavelet coherence; dynamic functional connectivity.

**Goals for Brainhack Marseille**
- Describe the methods to other participants
- Refine the method implementation to estimate spectral connectivity present in xfrites
- Create the documentation for the method
- Write smoke and functional unit tests
- Create examples illustrating the purpose of the single-trial coherence / PLV

**Skills:**
- Spectral analysis 50%
- Python 80%

**Striking Image**
![Dynamic functional connectivity estimates with single-trial coherence](image.png "coherence dFC")
<!-- Upload an image related to your project. -->

## Project submission

### Submission checklist

*Once the issue is submitted, please check items in this list as you add under 'Additional project info'*

- [x] Link to your project: could be a code repository, a shared document, etc.

- [x] Goals for the Brainhack: describe what you want to achieve during this brainhack. 

- [x] Skills: list skills that would be particularly suitable for your project (coding and non-coding)
![image](https://user-images.githubusercontent.com/17538901/142253801-e6a5ed54-e502-44db-84e7-6978fe993ef0.png)
.
