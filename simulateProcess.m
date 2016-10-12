% Function to simulate one choice drift diffusion model using a mixed
% simulation method as described in the paper:
%
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014. doi: 10.1145/2554850.2554872.
%
% RT = simulateProcess(theta,samples,maxRT)
%
% theta: parameter vector [a, Ter (sec), eta, St(sec), xi]
% a is the boundary, Ter is mean non-decision time, St is variability in
% non-decision time, xi is mean drift and eta is across trial variability
% in drift. Note that a/eta should have unit of seconds and eta/nu should be unitless. 
% typical value of theta is: [0.1,0.160,0.25,0.050,1.0]
% samples: is the number of simulated RTs required (default is 20,000)
% maxRT: (in seconds) is the maximum allowed RT in the PVT (default value is
% 30 seconds)
% (RT) returned are in seconds.
% (c)-2014 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

function RT = simulateProcess(theta,samples,maxRT)
%default parameter values

if (nargin<3)
    maxRT = 30; %the diffusion process must stop before maxRT
end   
if (nargin<2)
    samples = 20000; %more is better 
end

%s is within trial drift standard deviation
%Note that s is not a free parameter
s = 0.1;
tau=0.0005; %smaller is better

length = ceil(maxRT/tau);
delta=s*sqrt(tau);

%initialize some useful variables for speed
RT = zeros(samples,1);
index = false(samples,1);
mu = zeros(samples,1);
lambda = zeros(samples,1);
TerVec = zeros(samples,1);

for i=1:samples
    %drift is sampled from a normal distribution
    drift = normrnd(theta(5),theta(3));
    %Non decision time  is sampled from an uniform distribution
    Ter = theta(2)+(rand(1)-0.5)*theta(4);
    %If drift is > 0 then the process results in a shifted wald or inverse gaussian distribution
    %therefore we can easily sample from a wald distribution
    if(drift > 0)
        %store the parameters for wald distribution
        %the sampling will be done later in a single go 
        %this is due to speed considerations
        lambda(i,1) = (theta(1)/s)^2;
        mu(i,1) = theta(1)/drift;
        TerVec(i,1) = Ter;
        index(i) = true;
    else        
        %for negative drifts random walk approximation is the only way
        p=0.5*(1+drift*sqrt(tau)/s);
        %simple vectorization trick to carry out random walk in a single line 
        samplePath = [0 cumsum(2*delta.*(rand(1,length-1)<p)-delta)];
        hit = find(samplePath >= theta(1), 1, 'first');
        %check if there is no hit with boundary within maxRT
        if(isempty(hit))
            hit = length;
        end
        RT(i) = Ter + hit*tau;
    end
end

%By generating all samples for wald distribution in a single go the speed 
%of simulation increases by 20X or more.

RT(index) = shiftedWald(TerVec(index),mu(index),lambda(index));

%ceil the RT at max RT
RT(RT>maxRT)= maxRT;
end