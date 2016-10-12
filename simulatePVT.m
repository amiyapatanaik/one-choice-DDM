% MATLAB Code to simulate a psychomotor vigilance task (PVT)
%
% RT = simulatePVT(theta,time,maxRT)
%
% theta: parameter vector [a, Ter (sec), eta, St(sec), xi]
% a is the boundary, Ter is mean non-decision time, St is variability in
% non-decision time, xi is mean drift and eta is across trial variability
% in drift. Note that a/eta should have unit of seconds and eta/nu should be unitless. 
% typical value of theta is: [0.1,0.160,0.25,0.050,1.0]
% time is the length of PVT in minutes (default value is 10 minutes)
% maxRT (in seconds) is the maximum allowed RT in the PVT (default value is 10 seconds) 
% maxRT may have strong effect on the number of samples that may be generated within the specified time.
% Please note that time on task effects are NOT modeled here. reaction time
% (RT) returned are in seconds. Simulation done completely relying on random walks 
% as speed is not an issue here.
% Sample run: RT = simulatePVT([0.1,0.160,0.25,0.050,1.0]) 
%
% (c)-2014 Amiya Patanaik amiyain@gmail.com
%
% Please cite the following paper in any published results that rely on the code.
%
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014. doi: 10.1145/2554850.2554872.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.

function RT = simulatePVT(theta,time,maxRT)

if nargin < 3,
    maxRT = 10;
end
if nargin < 2,
    time = 10;
end


%within trial drift standard deviation
%Note that s is not a free parameter

s = 0.1;
tau=0.0005;
delta=s*sqrt(tau);
ln = ceil(maxRT/tau);

%infinite loop
i = 1;
totalTime = 0;

while 1
    startTime = rand(1)*8+2; %ISI
    drift = normrnd(theta(5),theta(3));
    Ter = theta(2)+(rand(1)-0.5)*theta(4);
    p=0.5*(1+drift*sqrt(tau)/s);
    samplePath = [0 cumsum(2*delta.*(rand(1,ln-1)<p)-delta)];
    hit = find(samplePath >= theta(1), 1, 'first');
    %check if there is no hit with boundary
    if(isempty(hit)) 
        hit = ln;
    end
    RT(i) = Ter + hit*tau;
    %ceil the RT at max RT
    if(RT(i) > maxRT),
        RT(i) = maxRT;
        %consider this to be a lapse
        totalTime = totalTime + startTime + maxRT;
        i = i+1;
        continue;
    end
    totalTime = totalTime + startTime + RT(i);
    if(totalTime > time*60) %Time over
        break;
    end
    i = i+1;
end

RT = RT';
end