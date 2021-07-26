# Simulated examples for identifying latent patient subgroups
The programs and data in this repository are used to identify latent subject subgroups using electronic health records (EHRs) data. 
Since we were not allowed to release our EHRs, we implemented proposed algorithm to simulated datasets for the purpose of illustration. 
Please refer to the following sections about how these files are used to illustrate our algorithm. 
The locations and summaries of the programs and data are provideds in this document. 
Also, you can find annotations and comments in R scripts.

## 1  Flow chart
![Web Figure 1](https://github.com/jitonglou/IdSubgroup/blob/main/README_figures/TRule_simdata_flowchart_1.png?raw=true "Web Figure 1")

## 2  Simulation and parameter estimation settings
In [this explanatory work](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_GenerateData_seed72.R), we simulated a dataset of two health markers for 10,000 subjects. 
For the *i*th subject, we assumed *Y<sub>i1</sub>(t)* was Gaussian distributed and *Y<sub>i2</sub>(t)* was Bernoulli distributed. 
Thus, the inverse functions of canonical link functions *g<sub>1</sub><sup>-1</sup>(z)=z* and *g<sub>2</sub><sup>-1</sup>(z)=e<sup>z</sup>/(1+e<sup>z</sup>)*. 
Since the distribution of *Y<sub>i1</sub>(t)* has a dispersion parameter, we set *&phi;<sub>i1</sub>(t)=0.1*.
We generated two covariates *X<sub>i1</sub> ~ Normal(0,1/3)* and *X<sub>i2</sub> ~ Bernoulli(0.5)-0.5*. 
Thus, *__X__<sub>i</sub>=(1, X<sub>i1</sub>, X<sub>i2</sub>)<sup>T</sup>* was a 3-dimensional vector of baseline covariates.</br> 
The maximum observation time *T<sub>i</sub>* for each subject was set to 12 (days). 
The measured time points for simulated markers were generated from two Poisson processes *dN<sub>ik</sub>(t)*, *k=1, 2*.
Their intensity functions are set to
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" title="\mathbb{E}\[dN_{i1}(t)|\boldsymbol{X}_i\]=\exp\{0.5X_{i1}+0.25X_{i2}+0.3L_{i11}(t)-0.1L_{i12}(t)\}dt" /></a>
</p>
&nbsp;&nbsp;&nbsp;&nbsp;and
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}&plus;0.25X_{i2}&plus;0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." title="\mathbb{E}\[dN_{i2}(t)|\boldsymbol{X}_i\]=1.2\exp\{0.5X_{i1}+0.25X_{i2}+0.3L_{i21}(t)-0.1L_{i22}(t)\}dt." /></a>
</p>

Thus, *__&gamma;__<sub>1</sub>=__&gamma;__<sub>2</sub>=(0.5, 0.25)<sup>T</sup>* and *__&eta;__<sub>1</sub>=__&eta;__<sub>2</sub>=(0.3, -0.1)<sup>T</sup>*.
Additionally, in the intensity functions, we let *L<sub>ik1</sub>(t)=1* if there exists measurements of *k*th marker in \[*t*-3,*t*); otherwise, *L<sub>ik1</sub>(t)=0*.
If *L<sub>ik1</sub>(t)=1*, then *L<sub>ik2</sub>(t)* is average value of all *Y<sub>ik</sub>(t)* in \[*t*-3,*t*); otherwise, *L<sub>ik2</sub>(t)=0*.</br>

The true values of *__&beta;__<sub>k</sub>(t)* were assumed to be
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation}\color{Black}{&space;\begin{pmatrix}&space;\boldsymbol{\beta}_1^T(t)&space;\\&space;\boldsymbol{\beta}_2^T(t)&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;-1.36&plus;\frac{t}{10}&space;&&space;\sin(0.76&plus;t)&space;&&space;\cos(-0.3&plus;t)&space;\\&space;\cos(-0.25&plus;t)&space;&&space;0.37&plus;\frac{t}{10}&space;&&space;\sin(-0.68&plus;t)&space;\end{pmatrix}.}&space;\end{equation*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation}\color{Black}{&space;\begin{pmatrix}&space;\boldsymbol{\beta}_1^T(t)&space;\\&space;\boldsymbol{\beta}_2^T(t)&space;\end{pmatrix}&space;=&space;\begin{pmatrix}&space;-1.36&plus;\frac{t}{10}&space;&&space;\sin(0.76&plus;t)&space;&&space;\cos(-0.3&plus;t)&space;\\&space;\cos(-0.25&plus;t)&space;&&space;0.37&plus;\frac{t}{10}&space;&&space;\sin(-0.68&plus;t)&space;\end{pmatrix}.}&space;\end{equation*}" title="\begin{equation}\color{Black}{ \begin{pmatrix} \boldsymbol{\beta}_1^T(t) \\ \boldsymbol{\beta}_2^T(t) \end{pmatrix} = \begin{pmatrix} -1.36+\frac{t}{10} & \sin(0.76+t) & \cos(-0.3+t) \\ \cos(-0.25+t) & 0.37+\frac{t}{10} & \sin(-0.68+t) \end{pmatrix}.} \end{equation*}" /></a>
</p>

Furthermore, we assumed the correlation structure of multivariate latent processes to be 
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\begin{equation*}\color{Black}{&space;\boldsymbol{\Omega}(t)=&space;\begin{pmatrix}&space;0.5&space;&&space;-0.25\\&space;-0.25&space;&&space;0.5&space;\end{pmatrix}&space;&plus;&space;\frac{1}{20}&space;\begin{pmatrix}&space;\sin(t&plus;2)&space;&&space;\cos(t-0.5)\\&space;\cos(t-0.5)&space;&&space;\sin(t&plus;3)&space;\end{pmatrix}}&space;\end{equation*}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?\begin{equation*}\color{Black}{&space;\boldsymbol{\Omega}(t)=&space;\begin{pmatrix}&space;0.5&space;&&space;-0.25\\&space;-0.25&space;&&space;0.5&space;\end{pmatrix}&space;&plus;&space;\frac{1}{20}&space;\begin{pmatrix}&space;\sin(t&plus;2)&space;&&space;\cos(t-0.5)\\&space;\cos(t-0.5)&space;&&space;\sin(t&plus;3)&space;\end{pmatrix}}&space;\end{equation*}" title="\begin{equation*}\color{Black}{ \boldsymbol{\Omega}(t)= \begin{pmatrix} 0.5 & -0.25\\ -0.25 & 0.5 \end{pmatrix} + \frac{1}{20} \begin{pmatrix} \sin(t+2) & \cos(t-0.5)\\ \cos(t-0.5) & \sin(t+3) \end{pmatrix}} \end{equation*}" /></a>
</p>
 &nbsp;&nbsp;&nbsp;&nbsp;and
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\operatorname{Cov}(\boldsymbol{\epsilon}(t),&space;\boldsymbol{\epsilon}(s))&space;=&space;\exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times&space;\frac{\boldsymbol{\Omega}(t)&plus;\boldsymbol{\Omega}(s)}{2}," target="_blank"><img src="https://latex.codecogs.com/gif.latex?\operatorname{Cov}(\boldsymbol{\epsilon}(t),&space;\boldsymbol{\epsilon}(s))&space;=&space;\exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times&space;\frac{\boldsymbol{\Omega}(t)&plus;\boldsymbol{\Omega}(s)}{2}," title="\operatorname{Cov}(\boldsymbol{\epsilon}(t), \boldsymbol{\epsilon}(s)) = \exp\left\{-\left\(\frac{t-s}{b}\right)^2\right\}\times \frac{\boldsymbol{\Omega}(t)+\boldsymbol{\Omega}(s)}{2}," /></a>
</p>

where *b=0.5*.</br>
In the [uploaded simulation dataset](https://github.com/jitonglou/IdSubgroup/blob/main/rdata/TRule_simdata_GenerateData_seed72_workspace.RData), some of the 10000 subjects has 0 measurement for both health markers. Thus, we excluded these subjects and the sample size for parameter estimation is 8338 subjects. We use the scaled Epanechnikov kernel
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{&plus;}" target="_blank"><img src="https://latex.codecogs.com/gif.latex?K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{&plus;}" title="K_{h_{1n}}(u)=\frac{3}{4h_{1n}}\left[1-\left\(\frac{u}{h_{1n}}\right\)^2\right\]_{+}" /></a>
</p>

as the kernel function to estimate *__&beta;__<sub>k</sub>(t)*. Similarly, the kernel function for estimating *__&Omega;__(t)* was set to the product of two scaled univariate Epanechnikov kernels
<p align="center">
<a href="https://www.codecogs.com/eqnedit.php?latex=\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{&plus;}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{&plus;}." target="_blank"><img src="https://latex.codecogs.com/gif.latex?\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{&plus;}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{&plus;}." title="\widetilde{K}_{h_{2n}}(u_1,u_2)=\frac{9}{16h_{2n}^2}\left\{1-\left\(\frac{u_1}{h_{2n}}\right\)^2\right\}_{+}\left\{1-\left\(\frac{u_2}{h_{2n}}\right\)^2\right\}_{+}." /></a>
</p>

We extended the data-adaptive method in [Cao et al. (2015)](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4643299/) and selected the optimal bandwidths among {0.1,0.2,...,0.5}. We found *h=0.3* for *__&beta;__<sub>k</sub>(t)* and *h=0.2* for *__&Omega;__(t)* were close to the optimal values of bandwidths. This set of bandwidth was used in the uploaded simulation data. Since the proposed estimation method was expected to have more stable performance at time points that not on two ends, we used the estimated *__&beta;__<sub>k</sub>(t)* and *__&Omega;__(t)* at time points *t=3,4,...,11* to identify latent subgroups.

## 3  Results
   - Following [TRule_simdata_GenerateData_seed72.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_GenerateData_seed72.R), one could obtain estimated *__&gamma;__<sub>k</sub>* and *__&eta;__<sub>k</sub>*, *k=1, 2*, as
   
   > | Marker | Parameter | True value | Estimator |
   > | ------ | --------- | ---------- | --------- |
   > | *Y<sub>1</sub>* | *&gamma;<sub>11</sub>* | 0.5  | 0.456  |
   > | Continuous      | *&gamma;<sub>12</sub>* | 0.25 | 0.229  |
   > |                 | *&eta;<sub>11</sub>*   | 0.3  | 0.297  |
   > |                 | *&eta;<sub>12</sub>*   | -0.1 | -0.102 |
   > | *Y<sub>2</sub>* | *&gamma;<sub>21</sub>* | 0.5  | 0.456  |
   > | Binary          | *&gamma;<sub>22</sub>* | 0.25 | 0.228  |
   > |                 | *&eta;<sub>21</sub>*   | 0.3  | 0.304  |
   > |                 | *&eta;<sub>22</sub>*   | -0.1 | -0.100 |

   - Following [TRule_simdata_est_case0_seed72_tp10.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_est_case0_seed72_tp10.R), one could obtain estimated *__&beta;__<sub>k</sub>(t)*, *k=1, 2*, and *__&Omega;__(t)* at *t=10* as
   > | Marker | Parameter | True value | Estimator |
   > | ------ | --------- | ---------- | --------- |
   > | *Y<sub>1</sub>* | *&beta;<sub>10</sub>*   | -0.360 | -0.363 |
   > | Continuous      | *&beta;<sub>11</sub>*   | -0.972 | -1.044 |
   > |                 | *&beta;<sub>12</sub>*   | -0.962 | -0.945 |
   > | *Y<sub>2</sub>* | *&beta;<sub>20</sub>*   | -0.948 | -0.923 |
   > | Binary          | *&beta;<sub>21</sub>*   | 1.370  | 1.418 |
   > |                 | *&beta;<sub>22</sub>*   | 0.105  | 0.093 |

 &nbsp;&nbsp;&nbsp;&nbsp;and
   
   > | Parameter | True value | Estimator |
   > | --------- | ---------- | --------- |
   > | *&sigma;<sub>11</sub>*  | 0.473 | 0.415 |
   > | *&sigma;<sub>22</sub>*  | 0.521 | 0.486 |
   > | *&sigma;<sub>12</sub>*  | -0.300 | -0.261 |

   - Following [TRule_simdata_DistMat_case0_seed72.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_DistMat_case0_seed72.R), one could obtain the hierarchical clustering results as Web Figure 3. We classified the 8338 subjects to 4 subgroups. However, subgroup 4 only has one subject, ID = 1814. We checked this subject and found it had an outlier measurement of health marker 1 at *t=6.744*.
   ![Web Figure 3](https://github.com/jitonglou/IdSubgroup/blob/main/README_figures/TRule_simdata_DistMat_rectdend.png?raw=true "Web Figure 3")

## 4  Programs
   - [codes/functions_simdata.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/functions_simdata.R)\
   R functions for generating simulated data.
   - [codes/functions_estparams.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/functions_estparams.R)\
   R functions for estimating parameters (intensity parameters *&gamma;<sub>k</sub>* and *&eta;<sub>k</sub>*, regression coefficients *&beta;<sub>k</sub>(t)*, covariance matrix of latent processes *&Omega;(t)*).
   - [codes/functions_DistMat.cpp](https://github.com/jitonglou/IdSubgroup/blob/main/codes/functions_DistMat.cpp)\
   C++ functions for calculating the similarity between each pair of subjects *S<sub>ij</sub>*.
   - [codes/TRule_simdata_GenerateData_seed72.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_GenerateData_seed72.R)\
   Explanatory R scripts for generating a simulated dataset and estimating intensity parameters. R random seed=72.
   - [codes/TRule_simdata_est_case0_seed72_tp10.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_est_case0_seed72_tp10.R)\
   Explanatory R scripts for estimating regression coefficients and the covariance matrix. Time point *t*=10, bandwidth *h<sub>1n</sub>=0.3*, bandwidth *h<sub>2n</sub>=0.2*.  
   - [codes/TRule_simdata_DistMat_case0_seed72.R](https://github.com/jitonglou/IdSubgroup/blob/main/codes/TRule_simdata_DistMat_case0_seed72.R)\
   Explanatory R scripts for calculating the similarity between subjects and identifying latent subgroups.
   
## 5  Data
   - [rdata/TRule_simdata_GenerateData_seed72_workspace.RData](https://github.com/jitonglou/IdSubgroup/blob/main/rdata/TRule_simdata_GenerateData_seed72_workspace.RData)\
   R workspace for estimating intensity parameters, regression coefficients, and covariance matrix. 
   - [rdata/TRule_simdata_GenerateData_seed72_intensity_est.RData](https://github.com/jitonglou/IdSubgroup/blob/main/rdata/TRule_simdata_GenerateData_seed72_intensity_est.RData)\
   R workspace of estimated intensity parameters *&gamma;<sub>k</sub>* and *&eta;<sub>k</sub>*.
   - [rdata/TRule_simdata_est_case0_seed72.RData](https://github.com/jitonglou/IdSubgroup/blob/main/rdata/TRule_simdata_est_case0_seed72.RData)\
   R workspace of estimated regression coefficients *&beta;<sub>k</sub>(t)* and covariance matrix *&Omega;(t)* across time points.
   - [rdata/TRule_simdata_DistMat_case0_seed72_dendrogram.RData](https://github.com/jitonglou/IdSubgroup/blob/main/rdata/TRule_simdata_DistMat_case0_seed72_dendrogram.RData)\
   R workspace of hierarchical clustering results based the similarity metric between each pair of subjects *S<sub>ij</sub>*.
   - [TRule_simdata_DistMat_rectdend.png](https://github.com/jitonglou/IdSubgroup/blob/main/README_figures/TRule_simdata_DistMat_rectdend.png)\
   Visualization of identified latent subgroups for 8338 subjects and 4 subgroups.