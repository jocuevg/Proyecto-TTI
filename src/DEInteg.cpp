//$Source$
//------------------------------------------------------------------------------
// DEInteg
//------------------------------------------------------------------------------
// Proyecto-TTI.
//
//
/**@file DEInteg.cpp
 * @brief Programacion de operacion DEInteg.
 *
 * @author Jose Cuevas Gil de Gomez
 * @bug No hay.
 */
//------------------------------------------------------------------------------
#include "..\include\DEInteg.hpp"
#include "..\include\sign_.hpp"

Matrix& DEInteg(Matrix& func (double t, Matrix& y), double t, double tout, double relerr, double abserr,int n_eqn, Matrix& y)
{
    // maxnum = 500;
    double twou  = 2*__DBL_EPSILON__;
    double fouru = 4*__DBL_EPSILON__;

    DE_STATE.DE_INIT = 1;       // Restart integration
    DE_STATE.DE_DONE = 2;       // Successful step
    DE_STATE.DE_BADACC = 3;     // Accuracy requirement could not be achieved
    DE_STATE.DE_NUMSTEPS = 4;   // Permitted number of steps exceeded
    DE_STATE.DE_STIFF = 5;      // Stiff problem suspected
    DE_STATE.DE_INVPARAM = 6;    // Invalid input parameters
        

    int State_ = DE_STATE.DE_INIT;
    bool PermitTOUT = true;         // Allow integration past tout by default
    double told = 0;

    // Powers of two (two(n)=2^n)
    Matrix& two  = zeros(14);
    two(1)=1.0;two(2)= 2.0;two(3)= 4.0;two(4)= 8.0;two(5)= 16.0;two(6)= 32.0;two(7)= 64.0; 
    two(8)=128.0;two(9)=256.0;two(10)= 512.0;two(11)= 1024.0;two(12)= 2048.0;two(13)= 4096.0;two(14)= 8192.0;

    Matrix& gstr = zeros(14); 
    gstr(1)=1.0;gstr(2)= 0.5;gstr(3)= 0.0833;gstr(4)= 0.0417;gstr(5)= 0.0264;gstr(6)= 0.0188;gstr(7)= 0.0143;
    gstr(8)= 0.0114;gstr(9)= 0.00936;gstr(10)= 0.00789;gstr(11)= 0.00679;gstr(12)= 0.00592;gstr(13)= 0.00524;gstr(14)= 0.00468;

    Matrix& yy    = zeros(n_eqn,1);    // Allocate vectors with proper dimension
    Matrix& wt    = zeros(n_eqn,1);
    Matrix& p     = zeros(n_eqn,1);
    Matrix& yp    = zeros(n_eqn,1);
    Matrix& phi   = zeros(n_eqn,17);
    Matrix& g     = zeros(14,1);
    Matrix& sig   = zeros(14,1);
    Matrix& rho   = zeros(14,1);
    Matrix& w     = zeros(13,1);
    Matrix& alpha = zeros(13,1);
    Matrix& beta  = zeros(13,1);
    Matrix& v     = zeros(13,1);
    Matrix& psi_  = zeros(13,1);

    // while(true)

    // Return, if output time equals input time

    if (t==tout)    // No integration
        return;
    

    // Test for improper parameters

    double epsilon = fmax(relerr,abserr);

    if ( ( relerr <  0.0         ) ||             // Negative relative error bound
        ( abserr <  0.0         ) ||             // Negative absolute error bound
        ( epsilon    <= 0.0         ) ||         // Both error bounds are non-positive
        ( State_  >  DE_STATE.DE_INVPARAM ) ||   // Invalid status flag
        ( (State_ != DE_STATE.DE_INIT) &&       
        (t != told)           ) ){
        State_ = DE_STATE.DE_INVPARAM;              // Set error code
        return;                                     // Exit
    }

    // On each call set interval of integration and counter for
    // number of steps. Adjust input error tolerances to define
    // weight vector for subroutine STEP.

    double del    = tout - t;
    double absdel = fabs(del);

    double tend   = t + 100.0*del;
    if (!PermitTOUT)
        tend = tout;
    

    int nostep = 0;
    int kle4   = 0;
    bool stiff  = false;
    double releps = relerr/epsilon;
    double abseps = abserr/epsilon;

    bool start;
    bool OldPermit=false;
    double x,h;
    double delsgn=0.0;

    if  ( (State_==DE_STATE.DE_INIT) || (!OldPermit) || (delsgn*del<=0.0) ){
        // On start and restart also set the work variables x and yy(*),
        // store the direction of integration and initialize the step size
        start  = true;
        x      = t;
        yy     = y;
        delsgn = sign_(1.0, del);
        h      = sign_( max(fouru*abs(x), abs(tout-x)), tout-x );
    }

    int kold;

    while (true){   // Start step loop

        // If already past output point, interpolate solution and return
        if (abs(x-t) >= absdel){
            Matrix& yout  = zeros(n_eqn,1);
            Matrix& ypout = zeros(n_eqn,1);
            g(2,1)   = 1.0;
            rho(2,1) = 1.0;
            double hi = tout - x;
            int ki = kold + 1;
            
            // Initialize w[*] for computing g[*]
            int temp1;
            for (int i=1;i<=ki;i++){
                temp1 = i;
                w(i+1,1) = 1.0/temp1;
            }
            // Compute g[*]
            double term = 0.0;
            double psijm1, gamma, eta;
            for (int j=2;j<=ki;j++){
                psijm1 = psi_(j,1);
                gamma = (hi + term)/psijm1;
                eta = hi/psijm1;
                for (int i=1;i<=ki+1-j;i++){
                    w(i+1,1) = gamma*w(i+1,1) - eta*w(i+2,1);
                }
                g(j+1,1) = w(2,1);
                rho(j+1,1) = gamma*rho(j,1);
                term = psijm1;
            }
            
            // Interpolate for the solution yout and for
            // the derivative of the solution ypout   
            int i;   
            for (int j=1;j<=ki;j++){
                i = ki+1-j;
                yout  = yout  + transpose(extract_column(phi,i+1)*g(i+1,1));
                ypout = ypout + transpose(extract_column(phi,i+1)*rho(i+1,1));
            }
            yout = y + yout*hi;
            y    = yout;
            State_    = DE_STATE.DE_DONE; // Set return code
            t         = tout;             // Set independent variable
            told      = t;                // Store independent variable
            OldPermit = PermitTOUT;
            return y;                       // Normal exit
        }                         
        
        // If cannot go past output point and sufficiently close,
        // extrapolate and return
        if ( !PermitTOUT && ( fabs(tout-x) < fouru*fabs(x) ) ){
            h = tout - x;
            yp = transpose(func(x,yy));          // Compute derivative yp(x)
            y = yy + yp*h;                // Extrapolate vector from x to tout
            State_    = DE_STATE.DE_DONE; // Set return code
            t         = tout;             // Set independent variable
            told      = t;                // Store independent variable
            OldPermit = PermitTOUT;
            return y;                       // Normal exit
        }
        
        // Test for too much work
        //   if (nostep >= maxnum)
        //       State_ = DE_STATE.DE_NUMSTEPS; // Too many steps
        //       if (stiff) 
        //           State_ = DE_STATE.DE_STIFF;// Stiffness suspected
        //       end
        //       y         = yy;                // Copy last step
        //       t         = x;
        //       told      = t;
        //       OldPermit = true;
        //       return;                        // Weak failure exit
        //   end
        
        // Limit step size, set weight vector and take a step
        h  = sign_(min(fabs(h), fabs(tend-x)), h);
        for (int l=1;l<=n_eqn;l++){
            wt(l,1) = releps*fabs(yy(l,1)) + abseps;
        }
        
        //   Step
        //                                                                   
        // Begin block 0                                                     
        //                                                                   
        // Check if step size or error tolerance is too small for machine    
        // precision.  If first step, initialize phi array and estimate a    
        // starting step size. If step size is too small, determine an       
        // acceptable one.                                                   
        //                                                                   

        if (fabs(h) < fouru*fabs(x)){
            h = sign_(fouru*fabs(x),h);
            bool crash = true;
            return y;           // Exit 
        }

        double p5eps  = 0.5*epsilon;
        bool crash  = false;
        g(2,1)   = 1.0;
        g(3,1)   = 0.5;
        sig(2,1) = 1.0;

        int ifail = 0;

        // If error tolerance is too small, increase it to an 
        // acceptable value.                                  

        double round = 0.0;
        for (int l=1;l<=n_eqn;l++){
            round = round + (y(l,1)*y(l,1))/(wt(l,1)*wt(l,1));
        }
        round = twou*sqrt(round);
        if (p5eps<round){
            epsilon = 2.0*round*(1.0+fouru);
            crash = true;
            return y;
        }

        int k;
        double hold;

        if (start){
            // Initialize. Compute appropriate step size for first step. 
            yp = func(x,y);
            double sum = 0.0;
            for (int l=1;l<=n_eqn;l++){
                phi(l,2) = yp(l);
                phi(l,3) = 0.0;
                sum = sum + (yp(l)*yp(l))/(wt(l,1)*wt(l,1));
            }
            sum  = sqrt(sum);
            double absh = fabs(h);
            if (epsilon<16.0*sum*h*h)
                absh=0.25*sqrt(epsilon/sum);
            
            h    = sign_(max(absh, fouru*abs(x)), h);
            hold = 0.0;
            double hnew = 0.0;
            k    = 1;
            kold = 0;
            start  = false;
            bool phase1 = true;
            bool nornd  = true;
            if (p5eps<=100.0*round){
                nornd = false;
                for (int l=1<=l;n_eqn;l++){
                    phi(l,16)=0.0;
                }
            }
        }

    //                                                                   
    // End block 0                                                       
    //                                                                   

    //                                                                   
    // Repeat blocks 1, 2 (and 3) until step is successful               
    //                                                                   
    while(true){
    
    //                                                                 
    // Begin block 1                                                   
    //                                                                 
    // Compute coefficients of formulas for this step. Avoid computing 
    // those quantities not changed when step size is not changed.     
    //                                                                 
    
    int kp1 = k+1;
    int kp2 = k+2;
    int km1 = k-1;
    int km2 = k-2;
    
    // ns is the number of steps taken with size h, including the 
    // current one. When k<ns, no coefficients change.           
    int ns;
    if (h !=hold)
        ns=0;
    
    if (ns<=kold)
        ns=ns+1;
    
    int nsp1 = ns+1;
    
    if (k>=ns){
        // Compute those components of alpha[*],beta[*],psi[*],sig[*] 
        // which are changed                                          
        beta(ns+1,1) = 1.0;
        int realns = ns;
        alpha(ns+1,1) = 1.0/realns;
        double temp1 = h*realns;
        sig(nsp1+1,1) = 1.0;
        if (k>=nsp1){
            for (int i=nsp1;i<=k;i++){
                int im1   = i-1;
                int temp2 = psi_(im1+1,1);
                psi_(im1+1,1) = temp1;
                beta(i+1,1)  = beta(im1+1,1)*psi_(im1+1,1)/temp2;
                temp1    = temp2 + h;
                alpha(i+1,1) = h/temp1;
                int reali = i;
                sig(i+2,1) = reali*alpha(i+1,1)*sig(i+1,1);
            }
        }
        psi_(k+1,1) = temp1;
        
        // Compute coefficients g[*]; initialize v[*] and set w[*].
        if (ns>1){
            // If order was raised, update diagonal part of v[*]
            if (k>kold){
                int temp4 = k*kp1;
                v(k+1,1) = 1.0/temp4;
                int nsm2 = ns-2;
                for (int j=1;j<=nsm2;j++){
                    int i = k-j;
                    v(i+1,1) = v(i+1,1) - alpha(j+2,1)*v(i+2,1);
                }
            }
            
            // Update V[*] and set W[*]
            int limit1 = kp1 - ns;
            double temp5  = alpha(ns+1,1);
            for (int iq=1;iq<=limit1;iq++){
                v(iq+1,1) = v(iq+1,1) - temp5*v(iq+2,1);
                w(iq+1,1) = v(iq+1,1);
            }
            g(nsp1+1,1) = w(2,1);
        }else{
            for (int iq=1;iq<=k;iq++){
                int temp3 = iq*(iq+1);
                v(iq+1,1) = 1.0/temp3;
                w(iq+1,1) = v(iq+1,1);
            }
        }
        
        // Compute the g[*] in the work vector w[*]
        int nsp2 = ns + 2;
        if (kp1>=nsp2){
            for (int i=nsp2;i<=kp1;i++){
                int limit2 = kp2 - i;
                double temp6  = alpha(i,1);
                for (int iq=1;iq<=limit2;iq++){
                    w(iq+1,1) = w(iq+1,1) - temp6*w(iq+2,1);
                }
                g(i+1,1) = w(2,1);
            }
        }
    } // if K>=NS
    
    //
    // End block 1
    //
    
    //
    // Begin block 2
    //
    // Predict a solution p[*], evaluate derivatives using predicted
    // solution, estimate local error at order k and errors at orders
    // k, k-1, k-2 as if constant step size were used.
    //   
    
    // Change phi to phi star
    if (k>=nsp1)
        for i=nsp1:k
            temp1 = beta(i+1);
            for l=1:n_eqn
                phi(l,i+1) = temp1 * phi(l,i+1);
            end
        end
    end
    
    // Predict solution and differences 
    for l=1:n_eqn
        phi(l,kp2+1) = phi(l,kp1+1);
        phi(l,kp1+1) = 0.0;
        p(l)       = 0.0;
    end
    for j=1:k
        i     = kp1 - j;
        ip1   = i+1;
        temp2 = g(i+1);
        for l=1:n_eqn
            p(l)     = p(l) + temp2*phi(l,i+1);
            phi(l,i+1) = phi(l,i+1) + phi(l,ip1+1);
        end
    end
    if (nornd)
        p = y + h*p;
    else
        for l=1:n_eqn
            tau = h*p(l) - phi(l,16);
            p(l) = y(l) + tau;
            phi(l,17) = (p(l) - y(l)) - tau;
        end
    end
    xold = x;
    x = x + h;
    absh = abs(h);
    yp = func(x,p);
    
    // Estimate errors at orders k, k-1, k-2 
    erkm2 = 0.0;
    erkm1 = 0.0;
    erk = 0.0;
    
    for l=1:n_eqn
        temp3 = 1.0/wt(l);
        temp4 = yp(l) - phi(l,1+1);
        if (km2> 0)
            erkm2 = erkm2 + ((phi(l,km1+1)+temp4)*temp3)...
                            *((phi(l,km1+1)+temp4)*temp3);
        end
        if (km2>=0)
            erkm1 = erkm1 + ((phi(l,k+1)+temp4)*temp3)...
                            *((phi(l,k+1)+temp4)*temp3);
        end
        erk = erk + (temp4*temp3)*(temp4*temp3);
    end
    
    if (km2> 0)
        erkm2 = absh*sig(km1+1)*gstr(km2+1)*sqrt(erkm2);
    end
    if (km2>=0)
        erkm1 = absh*sig(k+1)*gstr(km1+1)*sqrt(erkm1);
    end
    
    temp5 = absh*sqrt(erk);
    err = temp5*(g(k+1)-g(kp1+1));
    erk = temp5*sig(kp1+1)*gstr(k+1);
    knew = k;
    
    // Test if order should be lowered 
    if (km2 >0)
        if (max(erkm1,erkm2)<=erk)
            knew=km1;
        end
    end
    if (km2==0)
        if (erkm1<=0.5*erk)
            knew=km1;
        end
    end
    
    //
    // End block 2
    //
    
    //
    // If step is successful continue with block 4, otherwise repeat
    // blocks 1 and 2 after executing block 3
    //
    
    success = (err<=epsilon);
    
    if (~success)
    
        //
        // Begin block 3
        //
        
        // The step is unsuccessful. Restore x, phi[*,*], psi[*]. If
        // 3rd consecutive failure, set order to 1. If step fails more
        // than 3 times, consider an optimal step size. Double error
        // tolerance and return if estimated step size is too small
        // for machine precision.
        //
        
        // Restore x, phi[*,*] and psi[*]
        phase1 = false; 
        x = xold;
        for i=1:k
            temp1 = 1.0/beta(i+1);
            ip1 = i+1;
            for l=1:n_eqn
                phi(l,i+1)=temp1*(phi(l,i+1)-phi(l,ip1+1));
            end
        end
        
        if (k>=2)
            for i=2:k
                psi_(i) = psi_(i+1) - h;
            end
        end
        
        // On third failure, set order to one. 
        // Thereafter, use optimal step size   
        ifail = ifail+1;
        temp2 = 0.5;
        if (ifail>3) 
        if (p5eps < 0.25*erk)
            temp2 = sqrt(p5eps/erk);
        end
        end
        if (ifail>=3)
            knew = 1;
        end
        h = temp2*h;
        k = knew;
        if (abs(h)<fouru*abs(x))
            crash = true;
            h = sign_(fouru*abs(x), h);
            epsilon = epsilon*2.0;
            return;                 // Exit 
        end
        
        //
        // End block 3, return to start of block 1
        //
        
    end  // end if(success)
    
    if (success)
        break;
    end
    
    end

    //
    // Begin block 4
    //
    // The step is successful. Correct the predicted solution, evaluate
    // the derivatives using the corrected solution and update the
    // differences. Determine best order and step size for next step.
    //

    kold = k;
    hold = h;

    // Correct and evaluate
    temp1 = h*g(kp1+1);
    if (nornd)
        for l=1:n_eqn
            y(l) = p(l) + temp1*(yp(l) - phi(l,2));
        end
    else
        for l=1:n_eqn
            rho = temp1*(yp(l) - phi(l,2)) - phi(l,17);
            y(l) = p(l) + rho;
            phi(l,16) = (y(l) - p(l)) - rho;
        end
    end
    yp = func(x,y);

    // Update differences for next step 
    for l=1:n_eqn
        phi(l,kp1+1) = yp(l) - phi(l,2);
        phi(l,kp2+1) = phi(l,kp1+1) - phi(l,kp2+1);
    end
    for i=1:k
        for l=1:n_eqn
            phi(l,i+1) = phi(l,i+1) + phi(l,kp1+1);
        end
    end

    // Estimate error at order k+1 unless               
    // - in first phase when always raise order,        
    // - already decided to lower order,                
    // - step size not constant so estimate unreliable  
    erkp1 = 0.0;
    if ( (knew==km1) || (k==12) )
        phase1 = false;
    end

    if (phase1)
        k = kp1;
        erk = erkp1;
    else
        if (knew==km1)
            // lower order 
            k = km1;
            erk = erkm1;
        else
            if (kp1<=ns)
                for l=1:n_eqn
                    erkp1 = erkp1 + (phi(l,kp2+1)/wt(l))*(phi(l,kp2+1)/wt(l));
                end
                erkp1 = absh*gstr(kp1+1)*sqrt(erkp1);
                // Using estimated error at order k+1, determine 
                // appropriate order for next step               
                if (k>1)
                    if ( erkm1<=min(erk,erkp1))
                        // lower order
                        k=km1; erk=erkm1;
                    else
                        if ( (erkp1<erk) && (k~=12) )
                            // raise order 
                            k=kp1;
                            erk=erkp1;
                        end
                    end
                elseif (erkp1<0.5*erk)
                    // raise order
                    // Here erkp1 < erk < max(erkm1,ermk2) else    
                    // order would have been lowered in block 2.   
                    // Thus order is to be raised                  
                    k = kp1;
                    erk = erkp1;
                end
            end // end if kp1<=ns
        end // end if knew!=km1
    end // end if !phase1

    // With new order determine appropriate step size for next step
    if ( phase1 || (p5eps>=erk*two(k+2)) )
        hnew = 2.0*h;
    else
        if (p5eps<erk)
            temp2 = k+1;
            r = p5eps/erk^(1.0/temp2);
            hnew = absh*max(0.5, min(0.9,r));
            hnew = sign_(max(hnew, fouru*abs(x)), h);
        else
            hnew = h;
        end
    end
    h = hnew;

    //
    // End block 4
    //

    // Test for too small tolerances
    if (crash)
        State_    = DE_STATE.DE_BADACC;
        relerr    = epsilon*releps;       // Modify relative and absolute
        abserr    = epsilon*abseps;       // accuracy requirements
        y         = yy;                   // Copy last step
        t         = x;
        told      = t;
        OldPermit = true;
        return;                       // Weak failure exit
    end
    
    nostep = nostep+1;  // Count total number of steps
    
    // Count number of consecutive steps taken with the order of
    // the method being less or equal to four and test for stiffness
    kle4 = kle4+1;
    if (kold>  4)
        kle4 = 0;
    end
    if (kle4>=50)
        stiff = true;
    end
    
    end // End step loop
    
    //   if ( State_==DE_STATE.DE_INVPARAM )
    //       error ('invalid parameters in DEInteg');
    //       exit; 
    //   end
    //   if ( State_==DE_STATE.DE_BADACC )
    //       warning ('on','Accuracy requirement not achieved in DEInteg');
    //   end
    //   if ( State_==DE_STATE.DE_STIFF )
    //       warning ('on','Stiff problem suspected in DEInteg');
    //   end
    //   if ( State_ >= DE_STATE.DE_DONE )
    //       break;
    //   end
    //   
    // end
    return y;
    }