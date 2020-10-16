const mainlib = require("./main.js");

const opMatrix = (m1,m2,op)=> {
    
    if(! Array.isArray(m1[0])){
        m1 = m1.map((x) => [x]);
    }if(! Array.isArray(m2[0])){
        m2 = m2.map((x) => [x]);
    }
    if(m1.length !== m2.length || m1[0].length !==m2[0].length ){
        throw "Matrices aren't compatible"
    }
    const r = [];
    for (let i = 0; i < m1.length; i++) {
        r.push([]);
        const line = r[i];
        for (let j = 0; j < m1[0].length; j++) {
            switch(op){
                case 0:
                    line.push(m1[i][j] - m2[i][j]);    
                    break;
                case 1:
                    line.push(m1[i][j] + m2[i][j]);
                    break;
                case 2:
                    line.push(m1[i][j] * m2[i][j]);
                    break;
                case 3:
                    line.push(m1[i][j] /  m2[i][j]);
                    break;
                default:
                    throw "error";
                
            }
        }
    }
    //console.table(r);
    return r;
};
const addMatrix = (m1,m2) => opMatrix(m1,m2,1);
const subMatrix = (m1,m2)=> opMatrix(m1,m2,0);
const divMatrix = (m1,m2) => opMatrix(m1,m2,3);

const transposeMatrix = (m) => {
    const mt =[]; 
    for (let i = 0; i < m[0].length; i++) {
        mt.push(m.reduce((acc,cv)=> [...acc,cv[i]],[]));    
    }
    return mt;
};

;
//const result = bissectionMethod((x)=> x**2 - 4*Math.cos(x),0,10);
//console.log(result);

//const r = newtonMethod((x)=> x**2 - 4*Math.cos(x),
//(x)=>2*x+4*Math.sin(x),10);
//console.log(r);


//const r = newtonMethod((x)=> x**2 - 4,
//(x)=>2*x,10**-29);
//console.log(r);


//const result = newtonSecMethod((x)=> x**2 - 4*Math.cos(x),10);
//console.log(result);

//const r = inverseInterpolation((x)=> x**2 - 4*Math.cos(x),3,5,10);
//console.log(r);


//const r = systemNewtonMethod([(x)=> (x[0]+2*x[1]-2),
//(x) => (x[0]**2 + 4*x[1]**2 -4) ],
//[[(x)=> 1,(x) => 2],
//[(x) => 2*x[0],(x)=> (8*x[1])]],[2,3]);
//console.table(r);
//
//const r2 = systemBroydenMethod([(x)=> (x[0]+2*x[1]-2),
//(x) => (x[0]**2 + 4*x[1]**2 -4)],[[1,2],[4,24]],[2,3]);
//console.table(r2);

//const r = nonLinearLS((b,x,y)=> Math.exp(x**b[0]/b[1])-y
//,[(b,x)=> x**b[0]*(Math.log(x)/b[1])*Math.exp(x**b[0]/b[1]),
//(b,x) => -1*x**b[0]*(x**b[0]/b[1]**2)*Math.exp(x**b[0]/b[1])],[1,2,3],[1.995,1.410,1.260],[0,1]);
//console.log(r);


//const r =numericalIntegration((x)=>Math.exp(-(x**2)),0,1,10,false);
//console.log(r);

function bissectionMethod(f,a,b,tol=10**-3){
    if(!(f(a)<0 && f(b)>0)){
        throw "Incorrect Arguments";
    }
    let x,fv;
    do{
        x = (a+b)/2;
        fv = f(x);
        if(fv>0){
            b =x;
        }else{
            a =x;
        }
        
    }while(Math.abs(b-a) > tol);

    return x;
}

function newtonMethod(f,df,x,maxiter=100,tol=10**-3){
    for (let i=0;i<maxiter;i++) {
        let old_x = x;
        let derivative = df(old_x);
        if(derivative ===0 || derivative === Infinity){
            throw "Singularity Point";
        }
        x = old_x - f(old_x)/derivative;
        console.log(old_x,x,derivative);
        if(Math.abs(x-old_x)<tol){
            return x;
        }
    }
    throw "Convergence not reached";

}

function newtonSecMethod(f,x0,maxiter=100,tol=10**-3){
    let x1 = x0+10**-3;
    let fa = f(x0);
    for (let i=0;i<maxiter;i++) {
        let fi = f(x1);
        let i_derivative = ((x1-x0)/(fi-fa));
        x = x1 - fi*i_derivative;
        if(i_derivative ===0 || i_derivative === Infinity){
            throw "Singularity Point";
        }
        if(Math.abs(x-x1)<tol){
            return x;
        }
        [fa,x0,x1] = [fi,x1,x];
    }
    throw "Convergence not reached";

}

function inverseInterpolation(f,x0,x1,x2,maxiter=100,tol=10**-3){
    if(!(x0<x1<x2)){
        throw "Incorrect Parameters Order";
    }
    let old_x = Infinity;
    let [y1,y2,y3] = [f(x0),f(x1),f(x2)];
    for (let i = 0; i < maxiter; i++) {
        let x = (y2*y3*x0)/((y1-y2)*(y1-y3))
        + (y1*y3*x1)/((y2-y1)*(y2-y3))
        + (y1*y2*x2)/((y3-y1)*(y3-y2));
        if(Math.abs(x-old_x)<tol){
            return x;
        }
        old_x = x;
        switch (Math.max(y1,y2,y3)){
            case y1:
                x0 = old_x;
                [x0,x1,x2] = [x0,x1,x2].sort((a,b)=> a-b);
                [y1,y2,y3] = [f(x0),f(x1),f(x2)];
                break;
            case y2:
                x1 = old_x;
                [x0,x1,x2] = [x0,x1,x2].sort((a,b)=> a-b);
                [y1,y2,y3] = [f(x0),f(x1),f(x2)];
                break;
            case y3:
                x2 = old_x;
                [x0,x1,x2] = [x0,x1,x2].sort((a,b)=> a-b);
                [y1,y2,y3] = [f(x0),f(x1),f(x2)];
                break;
            default:
                throw "Unexpected Error"
        }
        
    }
    throw "Convergence not reached";
}


function systemNewtonMethod(fv,jmx,x0,maxiter=100,tol=10**-3){
    let x = x0;
    for (let i = 0; i < maxiter; i++) {
        let j,f,deltax;
        j = jmx.map((l)=>l.map((e) => e(x)));
        f = fv.map((f) => f(x));
        deltax = mainlib.solveSystem(j,f);
        x = x.map((ele,i)=> ele - deltax[i]);
        if((deltax.reduce((acc,x) => acc+x**2,0)/x.reduce((acc,x) => acc+x**2,0))**0.5 < tol){
            return x;
        }
    }
    throw "Convergence not reached";
}

function systemBroydenMethod(fv,b0,x0,maxiter=100,tol=10**-3){
    let x = x0;
    let b = b0;
    for (let i = 0; i < maxiter; i++) {
        let j,f,deltax,ex,oldFX,deltaxt;
        j = b;
        f = fv.map((f) => f(x));
        deltax = mainlib.solveSystem(j,f);
        x = x.map((ele,i)=> ele - deltax[i]);
        if((deltax.reduce((acc,x) => acc+x**2,0)/x.reduce((acc,x) => acc+x**2,0))**0.5 < tol){
            return x;
        }
        deltaxt = [deltax];
        ex = mainlib.matrixMultiplication(addMatrix(subMatrix(fv.map((f) => f(x)),f),
        mainlib.matrixMultiplication(b,deltax)),deltaxt);
        b = subMatrix(b,ex.map(l => l.map((e=> e/(mainlib.matrixMultiplication(deltaxt,deltax)[0][0])))));
    }
    throw "Convergence not reached";
}


function nonLinearLS(fcell,jcell,x,y,b0,maxiter=100,tol=10**-3){
    let b = b0;
    fv = x.map((v,i)=> ((b) => fcell(b,v,y[i])));
    console.log(fv);
    jmx = x.map((v)=> ((b) => jcell.map(f => f(b,v))));
    for (let i = 0; i < maxiter; i++) {
        let j,f,deltab,jtj,jt;
        j = jmx.map((g)=> g(b));
        f = fv.map((f) => f(b));
        jt = transposeMatrix(j);
        jtj = mainlib.matrixMultiplication(jt,j)
        deltab = mainlib.solveSystem(jtj,mainlib.matrixMultiplication(jt,f));
        b = b.map((ele,i)=> ele - deltab[i]);
        if((deltab.reduce((acc,x) => acc+x**2,0)/b.reduce((acc,x) => acc+x**2,0))**0.5 < tol){
            return b;
        }
    }
    throw "Convergence not reached";
}


function numericalIntegration(f,a,b,n,polyMethod=true){
    if( n >10 || n <1){
        throw "Unsupported numer of points"
    }
    if(polyMethod)
        return integrationPolynomialApproximation(f,a,b,n);
    return integrationGauss(f,a,b,n);

}

function integrationPolynomialApproximation(f,a,b,n){
    if(n<2)
        throw "Number of points is too small";
    let interval = (b-a)/(n-1);
    const points = [a];
    for (let i = 1; i < (n-1); i++) {
        points.push(a + interval*i);
    }
    points.push(b);
    const mvandermonde = [];
    for (let i = 0; i < n; i++) {
        mvandermonde.push([]);
        for (let j = 0; j < n; j++) {
            mvandermonde[i].push(points[j]**i);
        }
    }
    const bs = []
    for (let i = 1; i <= n; i++) {
        bs.push((b**i-a**i)/i)
    }
    const omega = mainlib.solveSystem(mvandermonde,bs);
    return omega.reduce((acc,cv,i)=>acc+cv*f(points[i]),0);
}

function integrationGauss(f,a,b,n){
    const gaussianQuadratureWeightsPoints = {
        1:[[2,0]],
        2:[[1.0000000000000000,-0.5773502691896257],[1.0000000000000000,0.5773502691896257]],
        3:[[0.8888888888888888,0.0000000000000000],[0.5555555555555556,-0.7745966692414834],[0.5555555555555556,0.7745966692414834]],
        4:[[0.6521451548625461,-0.3399810435848563],[0.6521451548625461,0.3399810435848563],[0.3478548451374538,-0.8611363115940526],[0.3478548451374538,0.8611363115940526]],
        5:[[0.5688888888888889,0.0000000000000000],[0.4786286704993665,-0.5384693101056831],[0.4786286704993665,0.5384693101056831],[0.2369268850561891,-0.9061798459386640],[0.2369268850561891,0.9061798459386640]],
        6:[[0.3607615730481386,0.6612093864662645],[0.3607615730481386,-0.6612093864662645],[0.4679139345726910,-0.2386191860831969],[0.4679139345726910,0.2386191860831969],[0.1713244923791704,-0.9324695142031521],[0.1713244923791704,0.9324695142031521]],
        7:[[0.4179591836734694,0.0000000000000000],[0.3818300505051189,0.4058451513773972],[0.3818300505051189,-0.4058451513773972],[0.2797053914892766,-0.7415311855993945],[0.2797053914892766,0.7415311855993945],[0.1294849661688697,-0.9491079123427585],[0.1294849661688697,0.9491079123427585]],
        8:[[0.3626837833783620,-0.1834346424956498],[0.3626837833783620,0.1834346424956498],[0.3137066458778873,-0.5255324099163290],[0.3137066458778873,0.5255324099163290],[0.2223810344533745,-0.7966664774136267],[0.2223810344533745,0.7966664774136267],[0.1012285362903763,-0.9602898564975363],[0.1012285362903763,0.9602898564975363]],
        9:[[0.3302393550012598,0.0000000000000000],[0.1806481606948574,-0.8360311073266358],[0.1806481606948574,0.8360311073266358],[0.0812743883615744,-0.9681602395076261],[0.0812743883615744,0.9681602395076261],[0.3123470770400029,-0.3242534234038089],[0.3123470770400029,0.3242534234038089],[0.2606106964029354,-0.6133714327005904],[0.2606106964029354,0.6133714327005904]],
        10:[[0.2955242247147529,-0.1488743389816312],[0.2955242247147529,0.1488743389816312],[0.2692667193099963,-0.4333953941292472],[0.2692667193099963,0.4333953941292472],[0.2190863625159820,-0.6794095682990244],[0.2190863625159820,0.6794095682990244],[0.1494513491505806,-0.8650633666889845],[0.1494513491505806,0.8650633666889845],[0.0666713443086881,-0.9739065285171717],[0.0666713443086881,0.9739065285171717]]
    }
    const l = (b-a);
    return l/2*gaussianQuadratureWeightsPoints[n]
    .map(([omega,z]) => [omega,(a+b+z*l)/2])
    .reduce((acc,[omega,x])=>acc+omega*f(x),0);

}

//console.log(numericalDerivative(Math.sin,1,"front"));
//console.log(numericalDerivative(Math.sin,1,"back"));
//console.log(derivativeRichard(Math.sin,1,3,"central",0.5));

function numericalDerivative(f,x,method,deltax){
    if(f(x)+10**-7 === f(x)){
        throw "deltax is too low"
    }
    switch(method){
        case "front":
            return (f(x+deltax)-f(x))/deltax;
        case "back":
            return (f(x)-f(x-deltax))/deltax;
        case "central":
            return (f(x+deltax)-f(x-deltax))/(2*deltax);
        default:
            throw "Unavailable method";
    }

}

function derivativeRichard(f,x,p,method,deltax=10**-2){
    const deltax2 = deltax/2;
    const q = deltax/deltax2;
    const d1 = numericalDerivative(f,x,method,deltax);
    const d2 = numericalDerivative(f,x,method,deltax2);
    return (d1 +(d1-d2)/(q**(-p)-1) )

}