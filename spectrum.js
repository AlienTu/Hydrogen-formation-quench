function qs(id){return document.getElementById(id)}
function range(a,b,n){return Array.from({length:n},(_,i)=>a+(b-a)*i/(n-1))}
function eig2(a,b,c){const tr=a+c,d=Math.sqrt((a-c)*(a-c)+4*b*b);return[(tr-d)/2,(tr+d)/2]}
function setStatus(s){qs('status').textContent=s}
function layout(title){return{title,margin:{t:36,l:45,r:15,b:35},paper_bgcolor:'#020617',plot_bgcolor:'#020617',font:{color:'#e5e7eb'}}}
function plotLine(id,traces,title){Plotly.newPlot(id,traces,layout(title),{responsive:true})}

function solveBackground(g,k1,k2,L,Nh){
  const sp=Math.sqrt(Math.PI), beta=2*sp, x=range(-L,0,Nh), dx=x[1]-x[0];
  const m1=Math.sqrt(2*Math.PI*k1), m2=Math.sqrt(2*Math.PI*k2);
  let p1=x.map(xx=>0.5*sp*(1+Math.tanh(0.5*m1*(xx+2/m1))));
  let p2=x.map(xx=>-0.5*sp*(1+Math.tanh(0.5*m2*(xx+2/m2))));
  p1[0]=0;p1[Nh-1]=sp/2;p2[0]=0;p2[Nh-1]=-sp/2;
  const dtau=0.06*dx*dx; let maxr=0;
  for(let it=0;it<900;it++){
    const n1=p1.slice(),n2=p2.slice();maxr=0;
    for(let i=1;i<Nh-1;i++){
      const d21=(p1[i-1]-2*p1[i]+p1[i+1])/(dx*dx), d22=(p2[i-1]-2*p2[i]+p2[i+1])/(dx*dx), ch=p1[i]+p2[i];
      const r1=-d21+0.5*g*g*ch+sp*k1*Math.sin(beta*p1[i]);
      const r2=-d22+0.5*g*g*ch+sp*k2*Math.sin(beta*p2[i]);
      n1[i]=p1[i]-dtau*r1; n2[i]=p2[i]-dtau*r2; maxr=Math.max(maxr,Math.abs(r1),Math.abs(r2));
    }
    p1=n1;p2=n2;p1[0]=0;p1[Nh-1]=sp/2;p2[0]=0;p2[Nh-1]=-sp/2;
  }
  return{x,p1,p2,dx,res:maxr}
}

function mirror(bg){const sp=Math.sqrt(Math.PI),h=bg.x,p1=bg.p1,p2=bg.p2;return{x:h.concat(h.slice(0,-1).map(v=>-v).reverse()),p1:p1.concat(p1.slice(0,-1).map(v=>sp-v).reverse()),p2:p2.concat(p2.slice(0,-1).map(v=>-sp-v).reverse()),dx:bg.dx,res:bg.res}}

function buildDiag(full,g,k1,k2){
  const beta=2*Math.sqrt(Math.PI),n=full.x.length-2, M11=new Float64Array(n), M22=new Float64Array(n), M12=0.5*g*g;
  for(let i=0;i<n;i++){const j=i+1;M11[i]=0.5*g*g+2*Math.PI*k1*Math.cos(beta*full.p1[j]);M22[i]=0.5*g*g+2*Math.PI*k2*Math.cos(beta*full.p2[j])}
  return{M11,M22,M12,n,dx:full.dx}
}

function applyH(v,out,op){
  const n=op.n,dx=op.dx,inv=1/(dx*dx),M12=op.M12,M11=op.M11,M22=op.M22;
  for(let i=0;i<n;i++){
    const u1=v[i],u2=v[i+n];
    const left1=i>0?v[i-1]:0,right1=i<n-1?v[i+1]:0,left2=i>0?v[i+n-1]:0,right2=i<n-1?v[i+n+1]:0;
    out[i]=(2*u1-left1-right1)*inv+M11[i]*u1+M12*u2;
    out[i+n]=(2*u2-left2-right2)*inv+M22[i]*u2+M12*u1;
  }
}
function dot(a,b){let s=0;for(let i=0;i<a.length;i++)s+=a[i]*b[i];return s}
function norm(a){return Math.sqrt(Math.max(dot(a,a),0))}
function axpy(y,a,x){for(let i=0;i<y.length;i++)y[i]+=a*x[i]}
function scale(y,a){for(let i=0;i<y.length;i++)y[i]*=a}

function lanczos(op,kdim){
  const dim=2*op.n, qPrev=new Float64Array(dim), q=new Float64Array(dim), z=new Float64Array(dim);
  for(let i=0;i<dim;i++) q[i]=Math.sin(12.9898*(i+1))+0.1*Math.cos(78.233*(i+1));
  scale(q,1/norm(q));
  const Q=[], alpha=[], beta=[];
  let b=0;
  for(let j=0;j<kdim;j++){
    Q.push(new Float64Array(q));
    applyH(q,z,op);
    if(j>0) axpy(z,-b,qPrev);
    const a=dot(q,z); alpha.push(a);
    axpy(z,-a,q);
    // one-step full reorthogonalization for stability
    for(let r=0;r<Q.length;r++){const c=dot(Q[r],z); axpy(z,-c,Q[r]);}
    b=norm(z); if(j<kdim-1) beta.push(b);
    if(b<1e-12) break;
    qPrev.set(q); q.set(z); scale(q,1/b);
  }
  return{Q,alpha,beta,m:alpha.length}
}

function tridiagEig(alpha,beta){
  const m=alpha.length,T=numeric.rep([m,m],0);
  for(let i=0;i<m;i++){T[i][i]=alpha[i]; if(i<m-1){T[i][i+1]=beta[i];T[i+1][i]=beta[i];}}
  const ev=numeric.eig(T), vals=ev.lambda.x, V=ev.E.x;
  const pairs=vals.map((v,i)=>({w2:v,w:Math.sqrt(Math.max(v,0)),i})).sort((a,b)=>a.w-b.w);
  return{pairs,V}
}

function reconstructVector(Q,coefCol){
  const dim=Q[0].length,u=new Float64Array(dim);
  for(let j=0;j<Q.length;j++){const c=coefCol[j];for(let i=0;i<dim;i++)u[i]+=c*Q[j][i];}
  const nm=norm(u); if(nm>0) scale(u,1/nm); return u;
}

function spectrumLanczos(full,g,k1,k2,num,kdim){
  const op=buildDiag(full,g,k1,k2); kdim=Math.min(kdim,2*op.n-2);
  const lan=lanczos(op,kdim), tri=tridiagEig(lan.alpha,lan.beta);
  const chosen=tri.pairs.slice(0,Math.min(num,tri.pairs.length));
  const n=op.n,dx=op.dx,iprs=[],w5=[],dens=[],vecs=[];
  for(const p of chosen){
    const coeff=lan.Q.map((_,j)=>tri.V[j][p.i]); const u=reconstructVector(lan.Q,coeff); vecs.push(u);
    let rho=[],nm=0; for(let r=0;r<n;r++){const rr=u[r]*u[r]+u[r+n]*u[r+n];rho.push(rr);nm+=rr*dx}
    rho=rho.map(v=>v/Math.max(nm,1e-300)); dens.push(rho);
    iprs.push(rho.reduce((s,v)=>s+v*v*dx,0)); w5.push(rho.reduce((s,v,r)=>s+(Math.abs(full.x[r+1])<5?v*dx:0),0));
  }
  return{pairs:chosen,iprs,w5,dens,vecs,krylov:lan.m}
}

function chooseModes(spec,lo,hi,topK){let chosen=[];spec.iprs.map((v,i)=>[v,i]).sort((a,b)=>b[0]-a[0]).slice(0,topK).forEach(p=>chosen.push(p[1]));const firstLo=spec.pairs.findIndex(m=>m.w>lo),firstHi=spec.pairs.findIndex(m=>m.w>hi);if(firstLo>=0)chosen.push(firstLo);if(firstHi>=0)chosen.push(firstHi);let upper=[];spec.pairs.forEach((m,i)=>{if(Math.abs(m.w-hi)<0.25)upper.push([spec.iprs[i],i])});upper.sort((a,b)=>b[0]-a[0]).slice(0,3).forEach(p=>chosen.push(p[1]));spec.pairs.forEach((m,i)=>{if(m.w<lo)chosen.push(i)});return Array.from(new Set(chosen)).sort((a,b)=>a-b)}

function drawBackgroundPlots(bg,g,k1,k2,lo2,hi2){
  const sp=Math.sqrt(Math.PI),pc=bg.p1.map((v,i)=>(v+bg.p2[i])/Math.sqrt(2)),ps=bg.p1.map((v,i)=>(v-bg.p2[i])/Math.sqrt(2));
  plotLine('backgroundPhi',[{x:bg.x,y:bg.p1,name:'phi1'},{x:bg.x,y:bg.p2,name:'phi2'}],'Background in phi1, phi2');
  plotLine('backgroundCS',[{x:bg.x,y:pc,name:'phi_c'},{x:bg.x,y:ps,name:'phi_s'}],'Background in phi_c, phi_s');
  const ml=[],mh=[];for(let i=0;i<bg.x.length;i++){const e=eig2(0.5*g*g+2*Math.PI*k1*Math.cos(2*sp*bg.p1[i]),0.5*g*g,0.5*g*g+2*Math.PI*k2*Math.cos(2*sp*bg.p2[i]));ml.push(e[0]);mh.push(e[1])}
  plotLine('massPlot',[{x:bg.x,y:ml,name:'lower local M2'},{x:bg.x,y:mh,name:'upper local M2'},{x:bg.x,y:bg.x.map(()=>lo2),name:'vac lower',line:{dash:'dot'}},{x:bg.x,y:bg.x.map(()=>hi2),name:'vac upper',line:{dash:'dot'}}],'Local mass-matrix eigenvalues')
}

function runSolver(){
  setStatus('Computing background...');
  setTimeout(()=>{try{
    const gin=+qs('g').value,k1=+qs('k1').value,k2=+qs('k2').value,L=+qs('L').value,Nh=Math.round(+qs('Nh').value),num=+qs('numModes').value,topK=+qs('topIpr').value,kdim=Math.round(+qs('krylovDim').value);
    const g=qs('sqrt2g').checked?Math.sqrt(2)*gin:gin;
    const vac=eig2(0.5*g*g+2*Math.PI*k1,0.5*g*g,0.5*g*g+2*Math.PI*k2),lo=Math.sqrt(vac[0]),hi=Math.sqrt(vac[1]);
    qs('lowerThr').textContent=lo.toFixed(4);qs('upperThr').textContent=hi.toFixed(4);
    const bg=mirror(solveBackground(g,k1,k2,L,Nh));qs('residual').textContent=bg.res.toExponential(2);drawBackgroundPlots(bg,g,k1,k2,vac[0],vac[1]);
    const dim=2*(bg.x.length-2); setStatus('Background plotted. Running matrix-free Lanczos: dim='+dim+', Krylov='+kdim+'...');
    setTimeout(()=>{try{
      const spec=spectrumLanczos(bg,g,k1,k2,num,kdim),sel=chooseModes(spec,lo,hi,topK);qs('selectedModes').textContent=sel.join(', ');
      Plotly.newPlot('spectrumPlot',[{x:spec.pairs.map((_,i)=>i),y:spec.pairs.map(m=>m.w),mode:'markers',marker:{size:spec.iprs.map(v=>7+260*v),color:spec.iprs,colorscale:'Viridis',showscale:true,colorbar:{title:'IPR'}}},{x:[0,spec.pairs.length-1],y:[lo,lo],mode:'lines',name:'lower threshold',line:{dash:'dot'}},{x:[0,spec.pairs.length-1],y:[hi,hi],mode:'lines',name:'upper threshold',line:{dash:'dot'}}],layout('Lanczos spectrum: marker size/color = IPR'),{responsive:true});
      const n=bg.x.length-2,xi=bg.x.slice(1,-1),modeTr=[],denTr=[];
      for(const j of sel){const m=spec.pairs[j],u=spec.vecs[j],arr=[];let mx=1e-12;for(let r=0;r<n;r++){const val=u[r]+u[r+n];arr.push(val);mx=Math.max(mx,Math.abs(val))}modeTr.push({x:xi,y:arr.map(v=>v/mx),name:'n='+j+', w='+m.w.toFixed(4)});const dm=Math.max(...spec.dens[j]);denTr.push({x:xi,y:spec.dens[j].map(v=>v/dm),name:'n='+j+', IPR='+spec.iprs[j].toFixed(3)})}
      plotLine('modePlot',modeTr,'Selected Lanczos-Ritz mode functions, all centered at 0');plotLine('densityPlot',denTr,'Selected mode densities');
      let html='<table><tr><th>n</th><th>omega</th><th>omega^2</th><th>IPR</th><th>W(|x|<5)</th><th>class</th></tr>';spec.pairs.forEach((m,i)=>{if(sel.includes(i)||i<8){const cl=m.w<lo?'below lower':(m.w<hi?'between':'above upper');html+=`<tr><td>${i}</td><td>${m.w.toFixed(5)}</td><td>${m.w2.toFixed(5)}</td><td>${spec.iprs[i].toFixed(4)}</td><td>${spec.w5[i].toFixed(3)}</td><td>${cl}</td></tr>`}});html+='</table>';qs('candidateTable').innerHTML=html;
      setStatus('Done. Lanczos used Krylov dimension '+spec.krylov+'. Increase Krylov dimension to improve high-index/upper-threshold candidates.');
    }catch(e){console.error(e);setStatus('Lanczos error: '+e.message)}} ,80)
  }catch(e){console.error(e);setStatus('Error: '+e.message)}},50)
}
qs('runBtn').onclick=runSolver;
