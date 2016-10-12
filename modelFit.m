% Computes a Chi2 fit between observed data and model
%
% For details refer to:
% Ratcliff, R., & Van Dongen, H. P. (2011). Diffusion model for one-choice reaction time
% tasks and the cognitive effects of sleep deprivation. Proceedings of the National
% Academy of Sciences, 108(27), 11285-11290.
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
function [chiStat,E,O] = modelFit(RTObs,theta,maxRT)
quantiles = 0:0.05:1;
iterations = 20000;
tau=0.0005;

if(nargin < 3),
    maxRT = max(RTObs);
end

RTObs(RTObs > maxRT) = maxRT;
bins = size(quantiles,2);
%Simulate process and obtain RTs
RT = simulateProcess(theta,iterations,maxRT);
observations = size(RTObs,1);
%Split data into bins
Quartiles = quantile(RTObs,quantiles');
chiStat = 0;
for i=2:bins,
    if(i==2),
         E(i-1) = (sum(RT < Quartiles(i))/size(RT,1))*observations;
         O(i-1) = (sum(RTObs < Quartiles(i))/size(RTObs,1))*observations;
    elseif(i==bins),
         E(i-1) = (sum(RT >= Quartiles(i-1))/size(RT,1))*observations;
         O(i-1) = (sum(RTObs >= Quartiles(i-1))/size(RTObs,1))*observations;
    else
         E(i-1) = (sum(RT < Quartiles(i) & RT >= Quartiles(i-1))/size(RT,1))*observations;
         O(i-1) = (sum(RTObs < Quartiles(i) & RTObs >= Quartiles(i-1))/size(RTObs,1))*observations;
    end
    
    if(E(i-1) == 0),
         chiStat = realmax;
       return;
    end
    chiStat = chiStat + ((O(i-1)-E(i-1))^2)/E(i-1);
end

