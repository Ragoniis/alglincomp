const multiplyAX = (matrix,x)=> matrix.reduce((acc,cv,i)=>{
    acc.push(cv.reduce((acc,cv,i)=> (acc+cv*x[i]),0));
    return acc;
},[]);
const idMatrix= (n)=>{
    const a = [];
    for (let i = 0; i < n; i++) {
        a.push([]);
        for (let j = 0; j < n; j++) {
            let value = (i===j)? 1:0;
            a[i].push(value);
        }
    }
    return a;
}
const matrixMultiplication = (a,b) =>{
    const c = [];
    let numColA,numColB;
    if(! Array.isArray(a[0])){
        a = a.map((x) => [x]);
    }if(! Array.isArray(b[0])){
        b = b.map((x) => [x]);
    }
    numColA = a[0].length;
    numColB =  b[0].length;
    if(!(numColA === b.length)){
        throw "Impossible to multiple this matrices";
    }
    for (let i=0;i<a.length;i++){
        let lineA = a[i];
        c.push([]);
        for (let j = 0; j < numColB; j++) {
            let sum =0;
            for(let k =0;k<b.length;k++){ 
                sum+= lineA[k]*b[k][j];
            }
            c[i].push(sum);
        }
    }
    return c;
}

matrix = [[5,-4,1,0],[-4,6,-4,1],[1,-4,6,-4],[0,1,-4,5]];
matrix2 = [[1,0.2,0],[0.2,1,0.5],[0,0.5,1]];
matrix3=[[3,-1,-1],[-1,3,-1],[-1,-1,3]];

matrix_c = [[16,9,8,7,6,5,4,3,2,1],[9,17,9,8,7,6,5,4,3,2],[8,9,18,9,8,7,6,5,4,3],
[7,8,9,19,9,8,7,6,5,4],[6,7,8,9,18,9,8,7,6,5],[5,6,7,8,9,17,9,8,7,6],[4,5,6,7,8,9,16,9,8,7],
[3,4,5,6,7,8,9,15,9,8],[2,3,4,5,6,7,8,9,14,9],[1,2,3,4,5,6,7,8,9,13]];

matrix_c2 = [[3,2,0],[2,3,-1],[0,-1,3]];

matrix_c3 = [[1,1,1,1],[1,2,4,8],[1,3,9,27],[1,4,16,64]];

matrix_p = [[22,4,6,10],[4,18,8,2],[6,8,28,10],[10,2,10,26]];
b_p = [5,5,6,6];

b_c3 = [1,2,9,20];

b_c = [4,0,8,0,12,0,8,0,4,0];
b = [-1,2,1,3];
b2 = [2,0,7,1];
b3=[1,2,1];
b_c2 = [1,-1,1];

//console.table(luDecomposition(matrix))
//solveSystem(matrix,b,true);
//console.log(determinant(matrix));
//cholesky(matrix);
//console.table(solveSystem(matrix_c3,b_c3));
//console.log(multiplyAX(matrix2,b2));
//console.log(solveSystem(matrix_c,b_c,false,true));
//console.log(determinant(matrix_c));
//console.log(powerMethod(matrix2));
//[m,x] =eignJacobi(matrix);
//console.table(x);
//console.table(m);
//console.table(matrixMultiplication(matrix2,matrix3));
//console.log(iterativeMethods(matrix,b,10**3));
//console.log(powerMethod(matrix_c2));
//[a,b] = eignJacobi(matrix_c2);
//console.table(a);
//console.table(b);
//console.log(iteratveMethods(matrix_c2,b_c2,10**-4,true));
//console.log(lineLinearRegression([1,2,3,4],[1,2.5,3.5,4.3]));
//console.table(cholesky(matrix_p));
//console.log(solveSystem(matrix_p,b_p,true));
//const xs = [1,2,3,4];
//const ys = [2,5,7,8.6];
//console.log(lineLinearRegression(xs,ys));
//[a,x] = eignJacobi(matrix_p);
//console.table(a);
//console.table(x);

function lineLinearRegression(xs,ys){
    const p = xs.map((cv)=>{
        return [1,cv];
    });
    const pt = [Array(xs.length).fill(1),xs.map((cv)=>{
        return cv;
    })];
    const a = matrixMultiplication(pt,p);
    const y = multiplyAX(pt,ys);
    const b = solveSystem(a,y);
    return b;

}

function eignJacobi(Matrix,tol=10**-8){
    x = Matrix.map((v,i) =>{
        return v.map((ov,oi)=>{
            return (i===oi) ? 1 : 0;
        });
    });
    Matrix = Matrix.slice();
    let count = 0;
    while(true){
        const max ={
            value:0,
            indexes:[],
            flag:false
        };
        for (let i = 0; i < Matrix.length; i++) {
            for(let j=0;j<Matrix.length;j++){
                let cv= Matrix[i][j];
                if((Math.abs(cv)>Math.abs(max.value)) && (i!==j)){
                    max.value = cv;
                    max.indexes= [i,j];
                    if(Matrix[i][i] == Matrix[j][j]){
                        max.flag = true;
                    }else{
                        max.flag = false;
                    }
                }
            }
        }
        if(Math.abs(max.value)<tol){
            break;
        }
        if(max.indexes[0] > max.indexes[1]){
            max.indexes = max.indexes.reverse();
        }
        let phi = (max.flag) 
        ? Math.PI/4 : 
        Math.atan((2*Matrix[max.indexes[0]][max.indexes[1]]/(Matrix[max.indexes[0]][max.indexes[0]]-Matrix[max.indexes[1]][max.indexes[1]])))/2;
        const p = idMatrix(Matrix.length);
        const pt= idMatrix(Matrix.length);

        p[max.indexes[0]][max.indexes[0]] = Math.cos(phi);
        p[max.indexes[0]][max.indexes[1]] = - Math.sin(phi);
        p[max.indexes[1]][max.indexes[0]] = Math.sin(phi);
        p[max.indexes[1]][max.indexes[1]] = Math.cos(phi);

        pt[max.indexes[0]][max.indexes[0]] = Math.cos(phi);
        pt[max.indexes[0]][max.indexes[1]] = Math.sin(phi);
        pt[max.indexes[1]][max.indexes[0]] = - Math.sin(phi);
        pt[max.indexes[1]][max.indexes[1]] = Math.cos(phi);
        x = matrixMultiplication(x,p);
        Matrix = matrixMultiplication(pt,matrixMultiplication(Matrix,p));
    }
    return [Matrix,x];
}

function iterativeMethods(Matrix,b,tol=10**-4,gs=false){
    for (let i = 0; i < Matrix.length; i++) {
        let di = Matrix[i][i];
        let sLine=0;
        let sCol=0;
        let j;
        for (let j = 0; j < Matrix.length; j++) {
            sLine+=Math.abs(Matrix[i][j]);
            sCol+=Math.abs(Matrix[j][i]);
            j=j;
        }
        sLine-=di;
        sCol-=di;
        if(sLine>di||sCol>di){
            throw "Matrix isn't diagonally dominant";
        }
    }
    const norm = (list)=>{
        return list.reduce((acc,cv)=>{
            return acc+cv**2;
        },0)**0.5;
    }
    let x = Matrix.map(()=> 1);
    let newx,old= [];
    let r = 1;
    do{
        if(gs){
            newx = x;
            old = x.slice();
        }
        else{
            newx = [];
        }
        for (let i = 0; i < Matrix.length; i++) {
            newx[i] = (b[i] - Matrix[i].reduce((acc,cv,index)=>{
                if(i===index)
                    return acc;
                return acc+cv*x[index];
            },0))/Matrix[i][i];
        }
        x = (gs)? old : x;
        r = (norm(newx.reduce((acc,cv,i)=> {
            acc.push(cv-x[i]);
            return acc;
        },[]))/norm(newx));
        x = newx;
    }while(r>=tol);
    return x;
}

function powerMethod(Matrix,tol=10**-3){
    let X = Matrix.map((_)=>1);
    let lambda,result;
    let newlambda=1;
    do{
        lambda = newlambda;
        result= multiplyAX(Matrix,X);
        newlambda = result[0];
        X = result.map((ele)=>(ele/newlambda));
    } while((Math.abs((newlambda-lambda)/(newlambda)))>=tol);
    return [newlambda,X];
}


function determinant(matrix){
    const lu = luDecomposition(matrix);
    return lu.reduce((acc,cv,i)=>{
        return acc*cv[i];
    },1);
}
function luDecomposition(A,doCopy=true){
    if(doCopy){
        matrix = A.map(l=>l.slice());
    }else{
        matrix= A;
    }
    size = matrix.length;
    for (p = 0; p <size; p++){
        if(matrix[p][p]===0){
            let changed = false;
            for (let line = p+1; line < matrix.length; line++) {
                if(matrix[line][p] !== 0){
                    changed = true;
                    temp = matrix[p];
                    matrix[p] =matrix[line];
                    matrix[line] =temp;
                    break 
                }
            }
            if(changed===false)
                throw "Determinant equals to 0.";
        }
        for(i = p+1;i<size;i++){
            matrix[i][p]/= matrix[p][p];
        }
        for(j = p+1;j<size;j++){
            for(k = p+1;k<size;k++){
                matrix[k][j] = matrix[k][j]-matrix[k][p]*matrix[p][j];  
            }
        }
    }
    return matrix;
}
function cholesky(matrix){
    const chol= matrix.map((el)=>[]);
    for (let index =0; index<matrix.length;index++){
        let sumSquareKnownLs = 0;
        for (let k = 0; k < index; k++) {
            sumSquareKnownLs+= chol[k][index]**2; 
        }
        const pdV = matrix[index][index]-sumSquareKnownLs;
        if(pdV < 0){
            throw "Matrix isn't positive definite";
        }
        chol[index][index] = (pdV)**0.5;
        for (let j = index+1; j < matrix.length; j++) {
            let sumKnownLs = 0;
            for (let k = 0; k < index; k++) {
                sumKnownLs+= chol[index][k]*chol[k][j]; 
            }
            chol[index][j] =(matrix[index][j]-sumKnownLs)/chol[index][index];
            chol[j][index] =chol[index][j];
            
        }
    }
    return chol;
}

function solveSystem(A,B,doCholesky=false,copyLu=true){
    const y = [];
    let lu;
    if(doCholesky){
        lu = cholesky(A);
        y[0] = B[0]/lu[0][0];
        for (let index = 1; index < lu.length; index++) {
            y[index] = (y.reduce((acc,cv,i)=>{
                return acc-cv*lu[index][i];
            },B[index]))/lu[index][index];
        }
    }
    else {
        lu = luDecomposition(A,copyLu);
        y[0] = B[0];
        for (let index = 1; index < lu.length; index++) {
            y[index] = y.reduce((acc,cv,i)=>{
                return acc-cv*lu[index][i];
            },B[index]);
        }
    }
    const x  = [];
    x[0] = y[lu.length-1]/lu[lu.length-1][lu.length-1];
    for (let index = (lu.length-2); index >= 0; index--) {
        x.push((x.reduce((acc,cv,i)=>{
            return acc-cv*lu[index][lu.length-1-i];
        },y[index]))/lu[index][index]);
    }
    x.reverse();
    return x;
}

module.exports = { solveSystem, matrixMultiplication };

function systemBroydenMethod(fv,b0,x0,maxiter=100,tol=10**-3){
    let x = x0;
    let b = b0;
    for (let i = 0; i < maxiter; i++) {
        let j,f,deltax,ex,oldFX,deltaxt;
        j = b;
        f = fv.map((f) => f(x));
        deltax = mainlib.solveSystem(j,f);
        oldFX = f(x);
        x = x.map((ele,i)=> ele - deltax[i]);
        if((deltax.reduce((acc,x) => acc+x**2,0)/x.reduce((acc,x) => acc+x**2,0))**0.5 < tol){
            return x;
        }
        deltaxt = [deltax];
        ex = mainlib.matrixMultiplication(subMatrix(subMatrix(f(x),oldFX),
        mainlib.matrixMultiplication(b,deltax)),deltaxt);
        b = addMatrix(b,ex.map(l => l.map((e=> e/(mainlib.matrixMultiplication(deltaxt,deltax)[0][0])))));
    }
    throw "Convergence not reached";
}