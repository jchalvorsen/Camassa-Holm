Content
=======

- Title: Exploring a finite difference scheme for the Camassa-Holm equation
- Abstract
- Introduction to Camassa-Holm
    - ~~Paper by Camassa and Holm, historical context~~
    - ~~Figures of solitons, peakons~~
    - ~~Multi-peakon solutions~~
    - ~~Peakon/antipeakon~~
    - ~~Attribute scheme to Holden/Raynaud~~
- Scheme
    - Figures, comparison to reference solution(s) (peakon)
    - CFL condition
    - Improvements as per diva-portal
    - Numerical convergence tests (space/time)
    - Errors
    - Justify correctness by evaluating behavior, compared to reference solution(s)
- Analysis
    - ~~Holden/Raynaud: main points of analysis~~
    - Consistency
    - Von Neumann: failure motivates linearization
    - "Matrix stability": unable to bound
- Linearization
    - Present scheme, and present modifications (and why?)
    - Motivation behind Taylor expansion about the peakon solution
    - Numerical results, compare with analytic
    - Analysis: Von neumann fails, "matrix stability" may be possible to prove, though too complex
- Implementation details
    - FFT/convolution for efficiency
    - memory/processing time "trade-off"
    - profiling?
- Conclusion
    - Thoughts on future experimentation/analysis
    - Elaborate on this once the rest of the report is semi-complete
 
