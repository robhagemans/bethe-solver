#include "roots.h"
#include "symmetric-string.h"
#include "solver.h"

using namespace std;

void showIteration (StateSolver* solver)
{
    cout<< "iter "<< solver->iterations() <<" ";
    cout<<" conv "<< solver->convergence() <<" ";
    cout<<" rts "<< getValues(solver->getRoots()) <<" ";
    cout<<endl;
};



int main()
{
    double desired_convergence = 1e-20;
    
    int N;
    StateSolver* all_roots = 0;
    
    cout<<"test: ";
    int test;
    cin>>test;
    
    cout<<" maxiter: ";
    int max_iter = 200;
    cin>>max_iter;
    
    //for(int test=1; test<=100;++test) 
    {
     
        cout<<endl<<"test case #"<<test<<endl;
           
        switch (test) {
        
        case 1:
            // only reals
            //iter 10  ( 0.216201306054027 0.598086994989464 ), i(  )  (  ), i(  )  conv 1.34717932379393e-23
            //#0#0: 2I [ +0 +2 -2 +4 -4 ;]  2J [+0 +2 -2 +4 -4] sums [+0;+2;-2;+4;-4;] sum +0 (+0,+0) (+0.21620131,+0) (-0.21620131,+0) (+0.59808699,+0) (-0.59808699,+0)
            N = 10;
            all_roots = new Bunch<SymRoots>({ new OriginRoot (), new RealPairs (N, {2, 4}) }, &showIteration);
            break;
            
        case 2:
                
            // 5-string
            //iter 23  (  ), i( 1.00001480038088 2.16376688302984 )  (  ), i(  )  conv 2.70076219416421e-24
            //#6#0: 2I [;;;; +0 ;]  2J [-8 -14 -2 -10 -6] sums [+0;] sum +0 (-3.1791694e-15,+2.1637669) (-1.8111943e-15,+1.0000148) (-1.8106887e-15,+0) (-1.8111943e-15,-1.0000148) (-3.1791694e-15,-2.1637669)
            N = 10;
            all_roots = new Bunch<SymRoots>( { new OriginRoot (), new ImagPairs (N, 0, {6, 8}) } , &showIteration);
            break;
            
        case 3:
            // 8-string
            // roots (0,0.5) (-0,-0.5) (0,1.50000218850042) (-0,-1.50000218850042) (0,2.52077526622762) (-0,-2.52077526622762) (0,4.49648402696482) (-0,-4.49648402696482)
            // qn (3.5,0) (4.5,0) (6.5,7.47447136404745e-12) (-6.5,-7.47465689640387e-12) (5.5,2.14258772906719e-11) (-5.5,-2.14259479696648e-11) (4.5,6.9150559699486e-14) (-4.5,-6.91682294477163e-14)
            N = 16;
            all_roots = new Bunch<SymRoots>( { new OriginPair (), new ImagPairs (N, 0, {9, 11, 13}) }, &showIteration );
            break;
        
        case 4:
      
            // single kite
            //#3#0: 2I [ +0 ;; +0 ;]  2J [-1 +9 +1 -9] sums [-1;+1;] sum +0 (+0.01414495,+0) (+0,+1.0035739) (-0.01414495,+0) (+0,-1.0035739)
            //iter 8  (  ), i(  )  ( 0.0141449501592471 ), i( 1.00357387588898 )  conv 2.27585135999056e-21
            N = 10;
            all_roots = new Bunch<SymRoots>(  { new Kite(N, 2) }, &showIteration );
            break;
            
        case 5:
            // kite + singular pair
            // 1,2,3
            //  roots [(0,0.5) (-0,-0.5) (0,1.58661641380926) (0.076669363041406,0) (-0,-1.58661641380926) (-0.076669363041406,0)]

            // qn [(2.5,0) (-2.5,0) (4.5,4.385e-10) (-0.5,0) (-4.5,-4.385e-10) (0.5,0)]
            // correct wth adjusted level
            {
                N = 12;
                Kite * c = new Kite(N, 2);
                c->setInitial(0.01, 0.501, 0.02, 0.502); 
                all_roots = new Bunch<SymRoots>( { new OriginPair(), c } , &showIteration );
                break;
            }
        case 6:
            // origin, singular pair, real pair
            //iter 5  ( 0.23612398353546 ), i(  )  (  ), i(  )  conv 1.76019828615352e-28
            //#1#0: 2I [ +0 +2 -2 ; +0 ;]  2J [+0 +2 -2 +4 +6] sums [+0;+2;-2;+10;] sum +10 (+4.2284582e-21,+0) (+0.23612398,+0) (-0.23612398,+0) (-2.1684043e-20,+0.5) (-2.1684043e-20,-0.5)
            N = 10;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new OriginPair(), new RealPairs (N, {2}) }, &showIteration );
            break;
        
        case 7:    
            // three string plus two real rapidities with zero qns
            //   roots (0,0) (0.119171104825923,0) (-0.119171104825923,0) (0,1.04460675205512) (-0,-1.04460675205512)
            //   qn (0,0) (3.23332361757864e-11,0) (-3.23332361757864e-11,0) (4,-6.01707173155083e-09) (-4,6.01707173155083e-09)
            //#3#0: 2I [ +1 -1 ;; +0 ;]  2J [+0 +0 +8 +2 +10] sums [+0;+0;+0;] sum +0 (+0.11917111,+0) (-0.11917111,+0) (+1.2598685e-16,+1.0446068) (+5.534098e-17,+0) (+1.2598685e-16,-1.0446068)
            //iter 100  ( 0.119171104825923 ), i( 1.04460675205512 )  (  ), i(  )  conv 2.56131024799912e-22
            N = 10;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new ImagPairs(N, 0, {8}), new RealPairs (N, {0}) } , &showIteration);
            break;
            
            
            
            
            
/** SymString test **/            
            
        case 8:
            // testing real roots through string logic

            N = 10;
            //#0#0: 2I [ +0 +2 -2 +4 -4 ;]  2J [+0 +2 -2 +4 -4] sums [+0;+2;-2;+4;-4;] sum +0 (+0,+0) (+0.21620131,+0) (-0.21620131,+0) (+0.59808699,+0) (-0.59808699,+0)
            all_roots = new Bunch<SymRoots>(  { new OriginRoot (), new SymString (N, 1, 2), new SymString (N, 1, 4) } , &showIteration );
            //    vector<Roots*> all_roots = { new OriginRoot (), new RealPairs (N, {2, 4}) } );
            break;
            
            
            
        case 9:
            // off-center 3-strings

            //iter 31  conv 4.48678e-21  rts [(0.676245,0.993633) (0.678017,-0) (0.676245,-0.993633) (-0.676245,0.993633) (-0.678017,0) (-0.676245,-0.993633)] 

            //roots [(0.676245041347245,0.99363339120859) (0.678017442184392,-0) (0.676245041347245,-0.99363339120859) (-0.676245041347245,0.99363339120859) (-0.678017442184392,0) (-0.676245041347245,-0.99363339120859)]

            //qn [(3.5,-6.56e-11) (2.5,0) (4.5,6.56e-11) (-4.5,-6.56e-11) (-2.5,0) (-3.5,6.56e-11)]

            N=12;
            // works, but what is the periodicity of 2xJ? this only works if I use 21, not e.g. 9 or -3...
            all_roots = new Bunch<SymRoots>( { new SymString(N, 3, 21) } , &showIteration );
            break;


         case 10:
            // 3&5-string kite

            //iter 39 lcnv 6.554e-19 scnv 7.784e-21 berr 6.487e-13
            //rts (0,2.352948977198935) (0.2267898966963307,1.00000212415697) (0.2267899397302771,0) (0.2267898966963307,-1.00000212415697) (-0,-2.352948977198935) (-0.2267898966963307,-1.00000212415697) (-0.2267899397302771,-0) (-0.2267898966963307,1.00000212415697)
            //qns  (7.5,5.695e-07) (5.5,1.095e-12) (-0.5,-1.414e-16) (6.5,-1.095e-12) (-7.5,-5.695e-07) (-5.5,-1.095e-12) (0.5,0) (-6.5,1.095e-12)

            N=16;
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {2.0}), new SymString(N, 3, 23) }, &showIteration );
//            all_roots = new Bunch<SymRoots>( { new SymString(N, 3, 23), new ImagPairs(N, {2.0})}, &showIteration );
  
            break;
            
         case 11:
            // 5&7-string kite

            N=24;
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {3.0}), new SymString(N, 5, 71-4*24) }, &showIteration );
        // roots (0,4.02628409383802) (0,-4.02628409383802) (0.526672708703533,2.00131375325592) (0.527668542662098,1.0000000002094) (0.527668542685181,0) (0.527668542662098,-1.0000000002094) (0.526672708703533,-2.00131375325592) (-0.526672708703533,2.00131375325592) (-0.527668542662098,1.0000000002094) (-0.527668542685181,-0) (-0.527668542662098,-1.0000000002094) (-0.526672708703533,-2.00131375325592) 

        // qn (11.5,7.47925103094376e-13) (-11.5,-7.48031121583758e-13) (9.49999999941701,-7.64913011080016e-10) (6.49999999055231,-5.87122328664888e-08) (1.50000002004641,3.53394964607057e-17) (7.49999999055232,5.87122326544518e-08) (10.499999999417,7.64913011080016e-10) (-10.499999999417,-7.64912640015304e-10) (-7.49999999055232,-5.87122332552233e-08) (-1.50000002004641,3.53394964607057e-17) (-6.49999999055232,5.87122331138653e-08) (-9.49999999941701,7.64912940401023e-10)

        // get jx2 from B-T i less the values on the imag axis! 
        // I thought: 11.5   9.5--10.5  7.5--8.5   -0.5  sum_j=35.5 -> 71
        // actual: 11.5  9.5--10.5  6.5--7.5  1.5 sum_j=35.5 -> 71

        // subtracting 2*24 from sum_jx2 does not work, but 4*24 does. is sum_j defined mod 2N?
            break;
        

        case 12:
            // 7&9-string kite

            N=32;
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {4.0}), new SymString(N, 7, 143) }, &showIteration );

            // I thought: 15.5   13.5--14.5  11.5--12.5  9.5--10.5 -0.5 sum_j = 71.5 -> 143 
            // actual: 15.5   13.5--14.5  8.5--9.5  8.5--9.5  7.5  sum_j = 71.5  
            // roots (0,5.95406368013793) (0,-5.95406368013793) (0.892620551124732,3.01724333719991) (0.899351665958059,1.99999940605389) (0.899349760222482,0.999999999999919) (0.899349760221857,0) (0.899349760222482,-0.999999999999919) (0.899351665958059,-1.99999940605389) (0.892620551124732,-3.01724333719991) (-0.892620551124732,3.01724333719991) (-0.899351665958059,1.99999940605389) (-0.899349760222482,0.999999999999919) (-0.899349760221857,-0) (-0.899349760222482,-0.999999999999919) (-0.899351665958059,-1.99999940605389) (-0.892620551124732,-3.01724333719991) 

            // qn (15.5,1.22161571365368e-12) (-15.5,-1.22152736491252e-12) (13.4999999999811,-3.56872724453755e-11) (8.49999999979562,-3.99492510580331e-11) (8.50000062794926,2.87895251814024e-06) (7.49999874453201,-1.2368823761247e-16) (9.50000062794927,-2.87895251879402e-06) (9.49999999979562,3.99488976630685e-11) (14.4999999999811,3.56873784638649e-11) (-14.4999999999811,-3.56869102155368e-11) (-9.49999999979562,-3.99486856260897e-11) (-9.50000062794927,2.87895251856431e-06) (-7.49999874453201,-1.2368823761247e-16) (-8.50000062794927,-2.8789525183346e-06) (-8.49999999979562,3.99486149470968e-11) (-13.4999999999811,3.56869543899074e-11)
            break;


        case 13:
// fails
            // 3&1-string kite
            {
            N=10;
            //all_roots = new Bunch<SymRoots>({ new ImagPairs(N, {1.0}), new SymString(N, 1, -1) } , &showIteration );
            // seems to converge, but extremely slowly.
            all_roots = new Bunch<SymRoots>({ new ImagPairs(N, {1.0}), new RealPairs(N, {-1}) } , &showIteration );
            // does not converge
            // the difference is in whether secant method (RealPairs) or simple iteration (SymString) is used
    
            // so far, the specialised Kite complex does much better.
            //all_roots = new Bunch<SymRoots>( { new Kite(N) } , &showIteration );
            break;
            }
            
        case 14:
            // 4&2-string kite

            cerr<<guessSum2xJ(12, { {}, {0}, {}, {0} }, 1, 0)<<endl; //8
            
// roots (0,1.63834422682211) (0,-1.63834422682211) (0.120258002930928,0.500000000016535) (0.120258002930928,-0.500000000016535) (-0.120258002930928,0.500000000016535) (-0.120258002930928,-0.500000000016535)

// qn (5.5,1.49878338439499e-12) (-5.5,-1.49874804489853e-12) (1.50000000000005,0.000260053029254424) (2.50000000000005,-0.000260053029254424) (-2.50000000000005,0.000260053029254282) (-1.50000000000005,-0.000260053029254282)
            N=12;
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {1.5}), new SymString(N, 2, 8) } , &showIteration );
            break;


        case 15:
            // 6&4-string kite
            
            //  roots (0,3.15528640172275) (0,-3.15528640172275) (0.364016507134246,1.50013308834115) (0.364126309841217,0.500000000000251) (0.364126309841217,-0.500000000000251) (0.364016507134246,-1.50013308834115) (-0.364016507134246,1.50013308834115) (-0.364126309841217,0.500000000000251) (-0.364126309841217,-0.500000000000251) (-0.364016507134246,-1.50013308834115)

            // qn (9.5,1.69412245108154e-12) (-9.5,-1.69428147881562e-12) (7.49999999992988,-7.47810249730879e-11) (2.50000000006558,1.17580980249108e-05) (3.50000000006558,-1.17580980248578e-05) (8.49999999992988,7.47811309915773e-11) (-8.49999999992988,-7.47810426428361e-11) (-3.50000000006558,1.17580980249991e-05) (-2.50000000006558,-1.17580980249285e-05) (-7.49999999992988,7.47812016705702e-11)


            N=20;
            cerr<<guessSum2xJ(20, { {}, {}, {}, {0}, {}, {0} }, 3, 0)<<endl; //44
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {2.5}), new SymString(N, 4, 44) } , &showIteration );
            // would expect 9.5    8.5--7.5  -0.5-- -0.5? -> 31
            break;

        case 16:
            // 3&7-string kite
            // needs 134 iterations
            // nimproved precision, now 239 :/ ... 
            // with more precise library, 111
            
            cerr<<guessSum2xJ(20, { {}, {}, {0}, {}, {}, {}, {0} }, 2, 0)<<endl; //25
            
            N = 20;
            all_roots = new Bunch<SymRoots>( { new ImagPairs(N, {2.0, 3.0}), new SymString(N, 3, 25) } , &showIteration );
        // 9.5 8.5    6.5--7.5  -1.5   --> 25
    //    roots (0,2.02613628926312) (0,3.81580028484744) (0,-3.81580028484744) (0,-2.02613628926312) (0.0567680255176861,0.999999997584167) (0.0567680217098974,0) (0.0567680255176861,-0.999999997584167) (-0.0567680255176861,0.999999997584167) (-0.0567680217098974,-0) (-0.0567680255176861,-0.999999997584167) 

    // qn (9.5,-2.80129191223715e-10) (8.5,1.04315125652711e-12) (-8.5,-1.04322193552003e-12) (-9.5,2.80129332581701e-10) (6.49999999670285,-2.02239414358073e-09) (-1.4999999934091,0) (7.49999999670285,2.02239428493871e-09) (-7.49999999670285,-2.02239446163619e-09) (1.4999999934091,0) (-6.49999999670285,2.02239456765468e-09)
            break;
        
        
        case 17:
            // off-center two-string
            
            // #0#0: 2I [; +1 -1 ;]  2J [+5 +5 -5 -5] sums [+10;+10;] sum +20 (+0.30544282,+0.49998451) (+0.30544282,-0.49998451) (-0.30544282,+0.49998451) (-0.30544282,-0.49998451)
            //  roots (0.305442823027925,0.499984508295454) (0.305442823027925,-0.499984508295454) (-0.305442823027925,0.499984508295454) (-0.305442823027925,-0.499984508295454)

            // qn (2.49999999995728,-3.07187960551211e-08) (2.49999999995728,3.07187960551211e-08) (-2.49999999995728,-3.07187959491026e-08) (-2.49999999995728,3.07187959491026e-08)

            N = 10;
            all_roots = new Bunch<SymRoots>( { new SymString(N, 2, 10) } , &showIteration );
       
            break;
            
        case 18:
            // off-center two-strings
            
            //#0#14: 2I [; +5 -5 ;]  2J [+7 +9 -9 -7] sums [-8;+8;] sum +0 (+1.2962217,+0.56352946) (+1.2962217,-0.56352946) (-1.2962217,+0.56352946) (-1.2962217,-0.56352946)
            // roots (1.29622170393265,0.563529463938434) (1.29622170393265,-0.563529463938434) (-1.29622170393265,0.563529463938434) (-1.29622170393265,-0.563529463938434)

            //qn (3.50000000000088,-3.14678955957014e-11) (4.50000000000088,3.14678955957014e-11) (-4.50000000000088,-3.14679132654496e-11) (-3.50000000000088,3.14679132654496e-11)
            
            N=12;
            all_roots = new Bunch<SymRoots>( { new SymString(N, 2, 16) } , &showIteration );

            break;
            
        case 19:
            // off-center two-strings

            // roots (0.686221122534119,0.498073185639976) (0.686221122534119,-0.498073185639976) (-0.686221122534119,0.498073185639976) (-0.686221122534119,-0.498073185639976)
            // qn (3.49999999996499,-1.22963648945864e-09) (3.49999999996499,1.22963648945864e-09) (-3.49999999996499,-1.22963650712838e-09) (-3.49999999996499,1.22963650712838e-09)

            N=12;
            all_roots = new Bunch<SymRoots>( { new SymString(N, 2, 14) } , &showIteration );

            break;




/** various kite attempts **/
        
        case 20:
            // 1,2,3 
// parity error
// not correct!

            N=12;
            all_roots = new Bunch<SymRoots>( { new OriginPair(), new Kite(N, 2) } , &showIteration );

            break;
        

        case 21:
            // 2,3 
// not correct! - parity error as imag pairs should be <1.0

            N=12;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new OriginPair(), new  ImagPairs(N, {1.0}), } , &showIteration );

            break;


        case 22:
            // 3, 4 
            // this is 1,6:
            // roots (0,0) (0,0.5) (0,-0.5) (0,2.86293271664476) (0,-2.86293271664476) (0,1.50089883934618) (0,-1.50089883934618)

            //qn (0,0) (3,0) (-3,0) (4,-4.71782277750422e-15) (-4,4.62947403635245e-15) (5,4.78249405602731e-13) (-5,-4.78337754343883e-13)


            N=14;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new OriginPair(), new ImagPairs(N, {1.0}), new ImagPairs(N, {1.5})} , &showIteration );

            break;
        case 23:
            // 5, 6 
            
//not correct!        
            
            //  roots (0,0) (0,0.5) (0,-0.5) (0,1.00000000000014) (0,2.00009988515701) (0,3.06177383575767) (0,5.5175373021582) (0,-5.5175373021582) (0,-3.06177383575767) (0,-2.00009988515701) (0,-1.00000000000014)

            //qn (0,5.65431943371292e-16) (5,0) (-5,0) (9.5,-0.000102743901965564) (8,-6.25302351300197e-11) (6.99999999999999,1.05877131396274e-13) (5.99999999999999,-1.45037710911794e-12) (-5.99999999999999,1.45034176962148e-12) (-6.99999999999999,-1.05912470892735e-13) (-8,6.25299877535444e-11) (-9.5,0.000102743901965069)
            
            // how do we differentiate 5,6 and 3,8 and 1,10 and 7,4 and 9,2?
            
            N=22;
            max_iter=400;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new OriginPair(), new ImagPairs(N, {1., 2., 3., 4.})} , &showIteration);

            break;

        case 24:
            // 1,2,4   --   4&2-string kite with origin root

            cerr<<guessSum2xJ(12, { {0}, {0}, {}, {0} }, 1, 0)<<endl; //8
            N=12;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new ImagPairs(N, {1.5}), new SymString(N, 2, 8) } , &showIteration );
            
            break;
            // roots (0,0) (0,2.20374560732863) (0,-2.20374560732863) (0.267110601848878,0.499999457969034) (0.267110601848878,-0.499999457969034) (-0.267110601848878,0.499999457969034) (-0.267110601848878,-0.499999457969034)

            //qn (0,0) (5,-6.69312393217536e-13) (-5,6.69224044476385e-13) (2.00000000000233,1.29092387690553e-08) (2.00000000000233,-1.29092384863393e-08) (-2.00000000000233,1.29092387337158e-08) (-2.00000000000233,-1.29092384509998e-08)
        
        
        case 25:
            // 1,3,5

            cerr<<guessSum2xJ(18, { {0}, {}, {0}, {}, {0} }, 2, 0)<<endl; //26
            
            // roots (0,0) (0,2.36702299920574) (0,-2.36702299920574) (0.235104646613796,1.00000191712805) (0.235106408083019,0) (0.235104646613796,-1.00000191712805) (-0.235104646613796,1.00000191712805) (-0.235106408083019,-0) (-0.235104646613796,-1.00000191712805)

            //qn (-7.06789929214115e-17,0) (8,8.82372450404308e-12) (-8,-8.82367149479839e-12) (5.99999989832656,-1.1061635966018e-07) (2.03323431466177e-07,0) (6.99999989832656,1.10616360084254e-07) (-6.99999989832656,-1.10616359907556e-07) (-2.03323431501517e-07,0) (-5.99999989832656,1.10616360296291e-07)


            N=18;
            all_roots = new Bunch<SymRoots>( { new OriginRoot(), new ImagPairs(N, {2.0}), new SymString(N, 3, 26) } , &showIteration );
            break;

        case 26:
            // also a 4-string, alternative using secant method (faster in this case)
            N=14;
            all_roots = new Bunch<SymRoots>( { new OriginPair(), new ImagPairs(N, {1.5}) } , &showIteration );
            break;
            
            
            
            
            
            
/** CentralString test **/            

        case 27:
            // 3, 4 
        
            // roots (0,1.00001564199356) (0,0) (-0,-1.00001564199356) (0,1.47803658780833) (0,0.5) (-0,-0.5) (-0,-1.47803658780833)

            //qn (5,-1.34122832135228e-08) (0,-3.53394964607057e-17) (-5,1.34122832135228e-08) (5,-1.86527375281053e-10) (3,0) (-3,0) (-5,1.86527375281053e-10)



            N=14;
            all_roots = new Bunch<SymRoots>( { new CentralString(N, 3), new CentralString(N, 4)} , &showIteration );

            break;
            
        case 28:
            // 4-string
            //  roots (0,1.5007358413045) (0,0.5) (-0,-0.5) (-0,-1.5007358413045)

            //qn (5.5,-1.05554552295884e-08) (3.5,0) (-3.5,0) (-5.5,1.05554553886161e-08)


            N=14;
            all_roots = new Bunch<SymRoots>( { new CentralString(N, 4)} , &showIteration );
            //all_roots = new Bunch<SymRoots>( { new OriginPair(), new ImagPairs({1.5}) } );

            break;
            
        case 29:
            // 5, 6 
            
            // TWO SOLUTIONS?
            
            // this one was found earlier with other sign rules...:
            // roots (0,1.99241303926923) (0,0.999999999987172) (0,0) (-0,-0.999999999987172) (-0,-1.99241303926923) (0,4.03026239417572) (0,1.50003837549915) (0,0.5) (-0,-0.5) (-0,-1.50003837549915) (-0,-4.03026239417572)

            //qn (8,3.13247794699096e-11) (9,1.55064156249267e-07) (0,-3.7106471283741e-16) (-9,-1.55064155948881e-07) (-8,-3.13249296627696e-11) (5.99999999999999,6.56445282556193e-12) (8,-2.15426036474816e-11) (5,0) (-5,0) (-8,2.15426389869781e-11) (-5.99999999999999,-6.56445282556193e-12)

            // current outcome:
            // roots (0,2.01697753200912) (0,1.00000000000872) (0,0) (-0,-1.00000000000872) (-0,-2.01697753200912) (0,2.28399476874115) (0,1.49998809737032) (0,0.5) (-0,-0.5) (-0,-1.49998809737032) (-0,-2.28399476874115)

            // qn (6.99999999999999,-4.43606070717679e-10) (9,9.15924211978471e-07) (0,1.76697482303529e-17) (-9,-9.1592421206682e-07) (-6.99999999999999,4.43605991203812e-10) (6.99999999999999,1.12900061179169e-11) (8,1.18826971818832e-10) (5,0) (-5,0) (-8,-1.188272898743e-10) (-6.99999999999999,-1.12898382553087e-11)

            N=22;
            all_roots = new Bunch<SymRoots>( { new CentralString(N, 5), new CentralString(N, 6)} , &showIteration );

            break;
        case 30:
            // 3, 8 
            // roots (0,1.00000000069283) (0,0) (-0,-1.00000000069283) (0,4.66145955588356) (0,2.54137773560629) (0,1.49999757531226) (0,0.5) (-0,-0.5) (-0,-1.49999757531226) (-0,-2.54137773560629) (-0,-4.66145955588356)

            // qn (9,-4.57814489316725e-09) (0,-6.18441188062351e-17) (-9,4.57814458394666e-09) (5.99999999999999,-1.44706579829957e-11) (6.99999999999999,1.10879702166511e-11) (9,9.18729724363082e-13) (5,0) (-5,0) (-9,-9.18420503769051e-13) (-6.99999999999999,-1.10879083725323e-11) (-5.99999999999999,1.44706756527439e-11)


            N=22;
            
            all_roots = new Bunch<SymRoots>( { new CentralString(N, 3), new CentralString(N, 8)} , &showIteration );

            break;

        case 31:
            // 8-string
            // roots (0,0.5) (-0,-0.5) (0,1.50000218850042) (-0,-1.50000218850042) (0,2.52077526622762) (-0,-2.52077526622762) (0,4.49648402696482) (-0,-4.49648402696482)
            // qn (3.5,0) (4.5,0) (6.5,7.47447136404745e-12) (-6.5,-7.47465689640387e-12) (5.5,2.14258772906719e-11) (-5.5,-2.14259479696648e-11) (4.5,6.9150559699486e-14) (-4.5,-6.91682294477163e-14)
            N = 16;
            
            all_roots = new Bunch<SymRoots>( { new CentralString(N, 8) } , &showIteration );
            break;
        
        case 32:
            // 7-string
            N = 14;
            // roots (0,3.62456588949357) (0,2.00539330578394) (0,1.00000000468672) (0,0) (-0,-1.00000000468672) (-0,-2.00539330578394) (-0,-3.62456588949357)

            //qn (4,-2.24376912487125e-11) (5,1.35385080848517e-11) (6,3.17777682034422e-10) (0,1.76697482303529e-17) (-6,-3.17777293299961e-10) (-5,-1.35385257545999e-11) (-4,2.24376470743419e-11)

            all_roots = new Bunch<SymRoots>( { new CentralString(N, 7) } , &showIteration );
            break;
        
        case 33:
            //5,8
            
            // roots (0,2.00086044097293) (0,1.00000000000001) (0,0) (-0,-1.00000000000001) (-0,-2.00086044097293) (0,5.11070502579276) (0,2.46047986198586) (0,1.49999990282132) (0,0.5) (-0,-0.5) (-0,-1.49999990282132) (-0,-2.46047986198586) (-0,-5.11070502579276)

            // qn (9,-5.35173383014224e-12) (11,0.00185177121515348) (0,-4.32908831643645e-16) (-11,-0.00185177121515279) (-9,5.3518840230022e-12) (6.99999999999999,-8.2916530453307e-12) (9,8.58898189880284e-12) (10,1.15599503955817e-10) (6,0) (-6,0) (-10,-1.15599062212111e-10) (-9,-8.5884871458524e-12) (-6.99999999999999,8.29170605457539e-12)

            N=26;
            desired_convergence = 1e-26;
            all_roots= new Bunch<SymRoots>( { new CentralString(N, 5), new CentralString(N, 8) } , &showIteration );
            break;
            

        case 34:
            // 2,3 
            
            // roots [(0,0.5) (0,-0.5) (0,0.999831) (0,0) (0,-0.999831)]

            // qn [(3,0) (-3,0) (5,-2.15218e-14) (0,0) (-5,2.14864e-14)]
 
            N=12;
            all_roots = new Bunch<SymRoots>( { new CentralString(N,2), new CentralString(N,3) } , &showIteration );

            break;            



/*** SymString coupling tests ***/


        case 44:
            // 3&5-string kite
            
            //iter 39 lcnv 6.554e-19 scnv 7.784e-21 berr 6.487e-13
            //rts (0,2.352948977198935) (0.2267898966963307,1.00000212415697) (0.2267899397302771,0) (0.2267898966963307,-1.00000212415697) (-0,-2.352948977198935) (-0.2267898966963307,-1.00000212415697) (-0.2267899397302771,-0) (-0.2267898966963307,1.00000212415697)
            //qns  (7.5,5.695e-07) (5.5,1.095e-12) (-0.5,-1.414e-16) (6.5,-1.095e-12) (-7.5,-5.695e-07) (-5.5,-1.095e-12) (0.5,0) (-6.5,1.095e-12)
            {
                N=16;
                
                CentralString* c = new CentralString(N, 5);
                SymString* d = new SymString(N, 3, 23);
                
                c->coupleString(*d, 1); 
                // simple iteration does best here
                c->setPolicy(0.4, 100);
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            }

        case 45:
            cerr<<guessSum2xJ(12, { {0}, {0}, {0} }, 0, 0)<<endl; //-1
            
            // 1&2&3-string kite
            // fails, see also case 43 (works) 20 (which fails).
            // WORKS with choice of initial values
            
            N=12;
            //roots [(0,0.5) (-0,-0.5) (0,1.58661641351762) (0,-1.58661641351762) (0.0766693630219695,0) (-0.0766693630219695,-0)]
            //qn [(2.5,0) (-2.5,0) (4.5,1.38e-12) (-4.5,-1.38e-12) (-0.5,0) (0.5,0)]

            {
                
                CentralString* c = new CentralString(N, 3);
                SymString* d = new SymString(N, 1, -1, 0.15);
                c->coupleString(*d, 2); 
                c->setPolicy(0.4,4);
                
                all_roots = new Bunch<SymRoots>( { new OriginPair(), c, d}, &showIteration );
                break;
            }
        case 50:
        
            // 7&9-string kite
// ok-ish, see also case 12, 40
            // negative argument to sqrt - I think this one is recoverable by taking sqrt(abs
            // however that introduces false positives
            // can probably be fixed by tuning initial values...
            N=32;
            {
                CentralString* c = new CentralString(N, 9);
                SymString* d = new SymString(N, 7, 143);
                c->coupleString(*d, 1); 
                c->setPolicy(1.5);
                
                all_roots = new Bunch<SymRoots>( { c, d  }, &showIteration );
                break;
            
            }
        case 51:
            // 5&7
            // ok, see also case 11, 41    
            N=24;
            {
                CentralString* c = new CentralString(N, 7);
                SymString* d = new SymString(N, 5, 71-4*24);
                c->coupleString(*d, 1); 
                // simple iteration does best here
                c->setPolicy(0.01,100);
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            
            }
        case 52:
            // 1&7
// fails (neg sqrt)
            // case 42 did not converge
            
            // roots [(0,4.04374625454729) (0,2.26369877071666) (0,1.25303416545882) (0,-1.25303416545882) (0,-2.26369877071666) (0,-4.04374625454729) (0.253034374906692,9.5367431640625e-27) (-0.253034374906692,-9.5367431640625e-27)]

            //qn [(5.5,-6.832e-11) (6.5,6.177e-10) (7.5,2.229) (-7.5,-2.229) (-6.5,-6.177e-10) (-5.5,6.832e-11) (-0.5,3.534e-17) (0.5,3.534e-17)]

            cerr<<guessSum2xJ(16, { {0}, {}, {}, {}, {}, {}, {0} }, 0, 0)<<endl;
            N=16;
            {
                CentralString* c = new CentralString(N, 7);
                SymString* d = new SymString(N, 1, -1, 1e-24);
                c->coupleString(*d, 1);
                c->setPolicy(0.1, 4);
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            
            }
        
        case 53:
            
            // 1&3-string kite
            // incorrect convergence
            // fails (neg sqrt)
            
            N=12;
           
            {
                
                CentralString* c = new CentralString(N, 3);
                SymString* d = new SymString(N, 1, -1,0.0001);
                c->coupleString(*d, 1); 
                
                c->setPolicy(0.1, 1);
                
                all_roots = new Bunch<SymRoots>( { c, d}, &showIteration );
                break;
            }
        
     case 54:    
            cerr<<guessSum2xJ(12, { {0}, {0}, {0} }, 0, 0)<<endl; //-1
            cerr<<guessSum2xJ(12, { {0}, {0}, {0} }, 2, 0)<<endl; //-1
            
            cerr<<guessSum2xJ(12, { {0}, {0}, {0} }, 0, 0)<<endl; //-1
            
            // 1&2&3-string kite
            // fails, see also case 43 (works) 20 (which fails).
            // WORKS  with choice of initial values
            // see case 45 for coprrect roots.
            
            N=12;
            
            {
                
                CentralString* c = new CentralString(N, 3);
                SymString* d = new SymString(N, 1, -1, 0.15); // initial epsilon, //-3
                c->coupleString(*d, 2); 
                c->setPolicy(0.4,4); //initial delta
                
                all_roots = new Bunch<SymRoots>( { new OriginPair(), c, d}, &showIteration );
                break;
            }
            
       case 55:
            cerr<<guessSum2xJ(20, { {}, {}, {0}, {}, {}, {0} }, 2, 0)<<endl; //28
            
            // 7&3-string kite
            // fails (neg sqrt)
            
            N=20;
            {
                
                CentralString* c = new CentralString(N, 7);
                SymString* d = new SymString(N, 3, 28);
                c->coupleString(*d, 1); 
                c->setPolicy(1.1);
                
                all_roots = new Bunch<SymRoots>( {  c, d}, &showIteration );
                break;
            }
             
            
            
            
            
/*** Kite coupling tests ***/            
            
       case 64:
            // 5&1-string kite N=16
            
            
            
             //roots [(0,2.01782280750237) (0,-2.01782280750237) (0,1.000000147823) (7.00194863471626e-05,0) (0,-1.000000147823) (-7.00194863471626e-05,0)]

            //qn [(6.5,6.306e-12) (-6.5,-6.306e-12) (7.5,1.751e-10) (-1.5,0) (-7.5,-1.751e-10) (1.5,0)]

            {
                N=16;
                desired_convergence = 1e-30;
                
                CentralString* c = new CentralString(N, 5);
                Kite* d = new Kite(N, 2);
                
                c->coupleKite(*d, 1); 
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            }
         case 65:
            // 5&1-string kite N=12
            //iter 18  conv 7.3827e-22  rts [(0,2.17773) (0,-2.17773) (0,1.00008) (0.00187936,0) (0,-1.00008) (-0.00187936,0)] 

            //roots [(0,2.17772809942668) (0,-2.17772809942668) (0,1.0000778692641) (0.00187935602790141,0) (0,-1.0000778692641) (-0.00187935602790141,0)]

            //qn [(4.5,-5.394e-12) (-4.5,5.394e-12) (5.5,-2.111e-11) (-1.5,0) (-5.5,2.111e-11) (1.5,0)]
            
            {
                N=12;
                
                CentralString* c = new CentralString(N, 5);
                Kite* d = new Kite(N, 2);
                
                c->coupleKite(*d, 1); 
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            }
        case 66:
            // 7&1-string kite
            // OK!
            
//            iter 28  conv 5.39726e-21  rts [(0,3.6378) (0,2.00592) (0,-2.00592) (0,-3.6378) (0,1) (-3.38443e-05,0) (0,-1) (3.38443e-05,0)] 

            //roots [(0,3.63779977695602) (0,2.00591511905848) (0,-2.00591511905848) (0,-3.63779977695602) (0,1.00000003473533) (-3.38443158324683e-05,0) (0,-1.00000003473533) (3.38443158324683e-05,0)]

            // qn [(5.5,-5.86e-12) (6.5,3.522e-12) (-6.5,-3.522e-12) (-5.5,5.859e-12) (7.5,-8.838e-12) (2.5,0) (-7.5,8.838e-12) (-2.5,0)]


            {
                N=16;
                
                CentralString* c = new CentralString(N, 7);
                Kite* d = new Kite(N, 2);
                
                c->coupleKite(*d, 1); 
                
                all_roots = new Bunch<SymRoots>( { c, d }, &showIteration );
                break;
            }
            
            
            
            
/** KiteString tests **/            
/*
        case 71:
            // 1&3-string kite
            // slow convergence (418 iter)
            
            N=10;
           {
                KiteString* c = new KiteString(N, 3);
                c->setPolicy(0.1);
                all_roots = new Bunch<SymRoots>( { c }, &showIteration );
                break;
            }
        case 72:
            
            // 3&5-string kite
            // no convergence
            
            N=16;
           {
                KiteString* c = new KiteString(N, 5);
                c->setPolicy(0.1);
                all_roots = new Bunch<SymRoots>( { c }, &showIteration );
                break;
            }
        case 73:
            
            // 1&2&3-string kite
            // parity error
          
            N=12;
           {
                KiteString* c = new KiteString(N, 3);
                c->setPolicy(0.502, 0.5);
                all_roots = new Bunch<SymRoots>( { new OriginPair(), c }, &showIteration );
                break;
            }
        
*/
        
/** String2 tests **/        

         case 80:
            // 3&5-string kite
// no convergence

            //rts (0,2.352948977198935) (0.2267898966963307,1.00000212415697) (0.2267899397302771,0) (0.2267898966963307,-1.00000212415697) (-0,-2.352948977198935) (-0.2267898966963307,-1.00000212415697) (-0.2267899397302771,-0) (-0.2267898966963307,1.00000212415697)
            //qns  (7.5,5.695e-07) (5.5,1.095e-12) (-0.5,-1.414e-16) (6.5,-1.095e-12) (-7.5,-5.695e-07) (-5.5,-1.095e-12) (0.5,0) (-6.5,1.095e-12)
            {
                N=16;
                SymString* c = new SymString(N, 5, -23);
                c->setOnAxis(1);
                c->setPolicy(1e-10);
                            
                all_roots = new Bunch<SymRoots>( { c }, &showIteration );
                break;
            }
         
         case 81:
            // 1&3-string kite
// incorrect convergence

            {
                N=10;
                SymString* c = new SymString(N, 3, -1);
                c->setOnAxis(1);
                c->setPolicy(1e-1);
                            
                all_roots = new Bunch<SymRoots>( { c }, &showIteration );
                break;
            }   
            
            
/** ImagPairs2 **/

            
        case 92:
            // 5-string

            // (  ), i( 1.00001480038088 2.16376688302984 )  (  ), i(  )  conv 2.70076219416421e-24
            //#6#0: 2I [;;;; +0 ;]  2J [-8 -14 -2 -10 -6] sums [+0;] sum +0 (-3.1791694e-15,+2.1637669) (-1.8111943e-15,+1.0000148) (-1.8106887e-15,+0) (-1.8111943e-15,-1.0000148) (-3.1791694e-15,-2.1637669)
            N = 10;
            all_roots = new Bunch<SymRoots>( { new OriginRoot (), new ImagPairs2 (N, {2, 4}) } , &showIteration);
            break;
            
        case 93:
            // 8-string

            // roots (0,0.5) (-0,-0.5) (0,1.50000218850042) (-0,-1.50000218850042) (0,2.52077526622762) (-0,-2.52077526622762) (0,4.49648402696482) (-0,-4.49648402696482)
            // qn (3.5,0) (4.5,0) (6.5,7.47447136404745e-12) (-6.5,-7.47465689640387e-12) (5.5,2.14258772906719e-11) (-5.5,-2.14259479696648e-11) (4.5,6.9150559699486e-14) (-4.5,-6.91682294477163e-14)
            N = 16;
            all_roots = new Bunch<SymRoots>( { new OriginPair (), new ImagPairs2 (N, {3, 5}) }, &showIteration );
            break;

        
/** Kite2 **/
        
        case 104:
      
            // single kite
            //#3#0: 2I [ +0 ;; +0 ;]  2J [-1 +9 +1 -9] sums [-1;+1;] sum +0 (+0.01414495,+0) (+0,+1.0035739) (-0.01414495,+0) (+0,-1.0035739)
            //iter 8  (  ), i(  )  ( 0.0141449501592471 ), i( 1.00357387588898 )  conv 2.27585135999056e-21
            N = 10;
            all_roots = new Bunch<SymRoots>(  { new Kite(N, 2) }, &showIteration );
            break;
            
        case 105:
            // kite + singular pair
            // 1,2,3
            //  roots [(0,0.5) (-0,-0.5) (0,1.58661641380926) (0.076669363041406,0) (-0,-1.58661641380926) (-0.076669363041406,0)]

            // qn [(2.5,0) (-2.5,0) (4.5,4.385e-10) (-0.5,0) (-4.5,-4.385e-10) (0.5,0)]
            // correct wth adjusted level
            {
                N = 12;
                Kite * c = new Kite(N, 2);
                c->setInitial(0.01, 0.501, 0.02, 0.502); 
                all_roots = new Bunch<SymRoots>( { new OriginPair(), c } , &showIteration );
                break;
            }
        case 106:
// no convergence
            // 3&5
            {
                N = 16;
                Kite* c = new Kite(N, 4);
                c->setInitial(0.2, 0.35, 0.2267899397302771, 0.352948977198935);
                SymString* d = new SymString(N, 3, 23);
                d->coupleKite(*c, 0);
                d->setPolicy(0.00000212415697, 0, 0.2);
                all_roots = new Bunch<SymRoots>( { c, d } , &showIteration );
                break;
            }

            //iter 39 lcnv 6.554e-19 scnv 7.784e-21 berr 6.487e-13
            //rts (0,2.352948977198935) (0.2267898966963307,1.00000212415697) (0.2267899397302771,0) (0.2267898966963307,-1.00000212415697) (-0,-2.352948977198935) (-0.2267898966963307,-1.00000212415697) (-0.2267899397302771,-0) (-0.2267898966963307,1.00000212415697)
            //qns  (7.5,5.695e-07) (5.5,1.095e-12) (-0.5,-1.414e-16) (6.5,-1.095e-12) (-7.5,-5.695e-07) (-5.5,-1.095e-12) (0.5,0) (-6.5,1.095e-12)

        }
        if (!all_roots) {
            cout<<" no test case with that number"<<endl;
            return 1;
        }
        
        IterResult solved = all_roots->solve(max_iter, desired_convergence);
        vector<CVal> roots = all_roots->getRoots();
        
        cout.precision(15);
        cout<<endl<<" test "<<test<<endl;
        cout<<" N "<<N<<endl; 
        cout<<" roots "<<getValues(roots)<<endl;
        cout.precision(4);
        cout<<endl<<" dirty "<< getDirtyJs(N, roots) <<endl;
        cout<<endl<<" clean "<< getJs(N, roots) <<endl;
        cout<<endl<<" BT "<< getJs_BT(N, roots) <<endl;

        double betherr = 0;
        bool parityerr = false;
        cout<<endl<<" jx2 "<< getJx2(N, roots, &betherr, &parityerr) <<endl;
        cout<<" err "<< betherr<<endl;
        cout<<endl;
        if (solved != IterResult::iter_ok) 
            cout<<"\033[1;31m "<<message(solved)<<" \033[0m"<<endl; 
        if (parityerr) 
            cout<<"\033[1;31m parity error \033[0m"<<endl; 
        if (betherr>1e-5)
            cout<<"\033[1;31m large quantum number error \033[0m"<<endl;
        
        delete all_roots;
        
    }

    return 0;    
}



