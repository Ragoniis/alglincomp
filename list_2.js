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


const r = systemNewtonMethod([(x)=> (x[0]+2*x[1]-2),
(x) => (x[0]**2 + 4*x[1]**2 -4) ],
[[(x)=> 1,(x) => 2],
[(x) => 2*x[0],(x)=> (8*x[1])]],[2,3]);
console.table(r);

const r2 = systemBroydenMethod([(x)=> (x[0]+2*x[1]-2),
(x) => (x[0]**2 + 4*x[1]**2 -4)],[[1,2],[4,24]],[2,3]);
console.table(r2);


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