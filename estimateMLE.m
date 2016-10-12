% Computes max likelihood estimate of DDM given the observed RTs (RTObs in
% seconds), tcensor is the censor time in seconds. The DDM parameters are not uniqely identifiable
% only boundary normalized drift ratios are unique. Renormalization is the boundary
% parameter for which you want the parameters estimated (or normalized) for.
%
% Tuned to work for standard PVTs. It is recommended that a mixed estimator
% be used for practical purposes. 
%
% For details refer to:
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014.
% doi: 10.1145/2554850.2554872.
%
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
function [theta, logl] = estimateMLE(RTObs,tcensor,renormalization)

    if (nargin<3)
        renormalization = 0.09; %boundary for normalization
    end 
	
	RTObs = sort(RTObs);
	eDrift = (renormalization*(0.1/std(RTObs))^2)^(1/3);
	eTer = mean(RTObs) - renormalization/eDrift;
	
	loss = @(p)loglikelihood(RTObs,tcensor,p');
	
	minVal = [0.05;0.11;0.01;0.01;0.2];
	maxVal = [0.12;0.30;0.8;0.15;2.0];
	
	initGuess = [renormalization,eTer,eDrift*0.5,eTer*0.25,eDrift];
	
    %check for bounds
	ubound = initGuess > maxVal';
	lbound = initGuess < minVal';
	
	initGuess(ubound) = maxVal(ubound')';
	initGuess(lbound) = minVal(lbound')';
	
	b = min(RTObs) - 0.001;
	a = [0,1,0,0.5,0];
	
	options = optimset('Algorithm','interior-point','TolCon',1e-15,'TolX',1e-15,'TolFun',1e-15,'FunValCheck','on','FinDiffType','central');
	[theta, logl] = fmincon(loss,initGuess,[],[],a,b,minVal,maxVal,[],options); 
	
	theta(5) = theta(5)*renormalization/theta(1);
	theta(3) = theta(3)*renormalization/theta(1);
	theta(1) = renormalization;

end