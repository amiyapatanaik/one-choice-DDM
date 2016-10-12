% Probablity density function of the one-choice drift diffusion model
% For details refer to:
% 
% Patanaik, A.,Zagorodnov, V. and Kwoh, C. K. Parameter estimation and simulation 
% for one-choice Ratcliff diffusion model., Symposium on Applied Computing-Association for Computing Machinery (SAC-ACM)
% Special Interest Group on Applied Computing (SIGAPP), 24-28 March 2014.
% doi: 10.1145/2554850.2554872.
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
function p = r(RT,theta)

    if(theta(4) == 0)
        p=0;
        return;
    end
        
    samples = size(RT,1);
    t1 = RT - repmat(theta(2)+theta(4)/2,samples,1);
    t2 = RT - repmat(theta(2)-theta(4)/2,samples,1);
    
    S1 = ones(samples,1);
    S2 = ones(samples,1);
    
    t1v = t1(t1>0);
    t2v = t2(t2>0);
      
    samples1 = size(t1v,1);
    samples2 = size(t2v,1);
    
    %s squared
    s = 0.01;
    %eta squared
    theta(3) = theta(3)*theta(3);  
      
    if(samples1 ~= 0)
    S1(t1>0) = normcdf((repmat(theta(1),samples1,1) - theta(5)*t1v)./realsqrt(theta(3)*(t1v.^2)+t1v*s),0,1)... 
        - exp((repmat(2*theta(1)*theta(5)+2*theta(1)*theta(1)*theta(3),samples1,1))./s)...
        .* normcdf(-(repmat(theta(1),samples1,1)*s + 2*theta(1)*theta(3)*t1v + theta(5)*t1v*s)./(realsqrt(theta(3)*t1v.^2+t1v*s)*s),0,1);
    end
    
    if(samples2 ~= 0)
    S2(t2>0) = normcdf((repmat(theta(1),samples2,1) - theta(5)*t2v)./realsqrt(theta(3)*t2v.^2+t2v*s),0,1)...
        - exp((repmat(2*theta(1)*theta(5)+2*theta(1)*theta(1)*theta(3),samples2,1))./s)...
        .* normcdf(-(repmat(theta(1),samples2,1)*s + 2*theta(1)*theta(3)*t2v + theta(5)*t2v*s)./(realsqrt(theta(3)*t2v.^2+t2v*s)*s),0,1);
    end
       
    p = (1/theta(4))*(S1-S2);
   
end