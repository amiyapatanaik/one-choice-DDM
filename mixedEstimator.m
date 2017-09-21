% Computes a mixed estimator of DDM given the observed RTs (RTObs in
% seconds), tcensor is the censor time in seconds. The DDM parameters are not uniqely identifiable
% only boundary normalized drift ratios are unique. Renormalization is the boundary
% parameter for which you want the parameters estimated (or normalized) for. This option is optional. 
% targetChiSq is the max ChiSq fit you are comfortable with, maxTries is
% maximum number of tries allowed.
%
% Estimator combines the standard MLE with Chi^2 based estimator
%
% It is recommended that this estimator is used for practical applications
%
% theta: parameter vector [a, Ter (sec), eta, St(sec), xi]
% Tuned to work for standard PVTs
%
% For details refer to:
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014.
% doi: 10.1145/2554850.2554872.
%
% AND
%
% Ratcliff, R., & Van Dongen, H. P. (2011). Diffusion model for one-choice reaction time
% tasks and the cognitive effects of sleep deprivation. Proceedings of the National
% Academy of Sciences, 108(27), 11285-11290.
%
% (c)-2014 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

function [theta, chiSq] = mixedEstimator(RTObs, tcensor, targetChiSq, maxTries, renormalization)

if (nargin<5)
    renormalization = 0.09; %boundary for normalization
end

if (nargin<3)
    targetChiSq = 30; %min ChiSq value required
end

if (nargin<4)
    maxTries = 4; %min ChiSq value required
end

%clear fast guesses
RTObs(RTObs<0.150) = []; 

initGuess = estimateMLE(RTObs,tcensor,renormalization);

minVal = [renormalization;0.11;0.01;0.01;0.2];
maxVal = [renormalization;0.30;0.8;0.15;2.0];

loss = @(p)modelFit(RTObs,p',tcensor);
chiSq = 100;
tries = 1;

while chiSq > targetChiSq && tries <= maxTries,
    fprintf('This is try number: %d',tries);
    [theta, chiSq] = fminsearchbnd(loss,initGuess,minVal,maxVal,optimset('TolFun',0.5,'TolX',0.05,'MaxFunEvals',150,'Display','iter'));
    initGuess = theta;
    tries = tries + 1;
end
end
