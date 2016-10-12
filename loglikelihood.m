% Log likelihood of data, for details refer to:
% 
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014.
% doi: 10.1145/2554850.2554872.
%
% RT: vector of observed reaction time in seconds
% theta: parameter vector [a, Ter (sec), eta, St(sec), xi]
% a is the boundary, Ter is mean non-decision time, St is variability in
% non-decision time, xi is mean drift and eta is across trial variability
% in drift. Note that a/eta should have unit of seconds and eta/nu should be unitless. 
% typical value of theta is: [0.1,0.160,0.25,0.050,1.0]
% tcensor: (in seconds) is the censor time  
% (c)-2014 Amiya Patanaik amiyain@gmail.com
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
% THE SOFTWARE.
function logl = loglikelihood(RT,tcensor,theta)
   count =  sum(RT>=tcensor);
   R = ddmCDF(tcensor,theta);
   RT(RT>=tcensor) = [];
   logl = -sum(log(r(RT,theta)))-count*log(1-R);
   if(isinf(logl) || isnan(logl) || ~isreal(logl)),
       logl = realmin;
   end
end