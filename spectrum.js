function qs(id){return document.getElementById(id)}
function range(a,b,n){return Array.from({length:n},(_,i)=>a+(b-a)*i/(n-1))}
function eig2(a,b,c){const tr=a+c, d=Math.sqrt((a-c)*(a-c)+4*b*b);return [(tr-d)/2,(tr+d)/2]}
function setStatus(s){qs('status').textContent=s}

function solveBackground(g,k1,k2,L,Nh){
  const sp=Math.sqrt(Math.PI), beta=2*sp, x=range(-L,0,Nh), dx=x[1]-x[0];
  const m1=Math.sqrt(2*Math.PI*k1), m2=Math.sqrt(2*Math.PI*k2);
  let p1=x.map(xx=>0.5*sp*(1+Math.tanh(0.5*m1*(xx+2/m1))));
  let p2=x.map(xx=>-0.5*sp*(1+Math.tanh(0.5*m2*(xx+2/m2))));
  p1[0]=0; p1[Nh-1]=sp/2; p2[0]=0; p2[Nh-1]=-sp/2;
  const dtau=0.10*dx*dx;
  let maxr=0;
  for(let it=0;it<900;it++){
    const n1=p1.slice(), n2=p2.slice(); maxr=0;
    for(let i=1;i<Nh-1;i++){
      const d21=(p1[i-1]-2*p1[i]+p1[i+1])/(dx*dx);
      const d22=(p2[i-1]-2*p2[i]+p2[i+1])/(dx*dx);
      const ch=p1[i]+p2[i];
      const r1=-d21+0.5*g*g*ch+sp*k1*Math.sin(beta*p1[i]);
      const r2=-d22+0.5*g*g*ch+sp*k2*Math.sin(beta*p2[i]);
      n1[i]=p1[i]-dtau*r1; n2[i]=p2[i]-dtau*r2;
      maxr=Math.max(maxr,Math.abs(r1),Math.abs(r2));
    }
    p1=n1; p2=n2; p1[0]=0; p1[Nh-1]=sp/2; p2[0]=0; p2[Nh-1]=-sp/2;
  }
  return {x,p1,p2,dx,res:maxr};
}

function mirror(bg){
  const sp=Math.sqrt(Math.PI), h=bg.x, p1=bg.p1, p2=bg.p2;
  const xr=h.concat(h.slice(0,-1).map(v=>-v).reverse());
  const f1=p1.concat(p1.slice(0,-1).map(v=>sp-v).reverse());
  const f2=p2.concat(p2.slice(0,-1).map(v=>-sp-v).reverse());
  return {x:xr,p1:f1,p2:f2,dx:bg.dx,res:bg.res};
}

function hessian(full,g,k1,k2){
  const beta=2*Math.sqrt(Math.PI), x=full.x, n=x.length-2, dx=full.dx;
  const H=numeric.rep([2*n,2*n],0), inv=1/(dx*dx);
  for(let i=0;i<n;i++){
    const j=i+1;
    const M11=0.5*g*g+2*Math.PI*k1*Math.cos(beta*full.p1[j]);
    const M22=0.5*g*g+2*Math.PI*k2*Math.cos(beta*full.p2[j]);
    H[i][i]=2*inv+M11; H[i+n][i+n]=2*inv+M22; H[i][i+n]=0.5*g*g; H[i+n][i]=0.5*g*g;
    if(i>0){H[i][i-1]=-inv; H[i+n][i+n-1]=-inv}
    if(i<n-1){H[i][i+1]=-inv; H[i+n][i+n+1]=-inv}
  }
  return H;
}

function spectrum(full,g,k1,k2,num){
  const H=hessian(full,g,k1,k2), ev=numeric.eig(H), vals=ev.lambda.x, V=ev.E.x;
  const pairs=vals.map((v,i)=>({w2:v,w:Math.sqrt(Math.max(v,0)),i})).sort((a,b)=>a.w-b.w).slice(0,num);
  const n=full.x.length-2, dx=full.dx, ipr=[], w5=[], dens=[];
  for(const m of pairs){
    let rho=[], norm=0;
    for(let r=0;r<n;r++){const u1=V[r][m.i], u2=V[r+n][m.i], rr=u1*u1+u2*u2; rho.push(rr); norm+=rr*dx}
    rho=rho.map(v=>v/norm); dens.push(rho);
    ipr.push(rho.reduce((s,v)=>s+v*v*dx,0));
    w5.push(rho.reduce((s,v,r)=>s+(Math.abs(full.x[r+1])<5?v*dx:0),0));
  }
  return {pairs,V,ipr,w5,dens};
}

function chooseModes(spec,lo,hi,topK){
  let chosen=[];
  spec.ipr.map((v,i)=>[v,i]).sort((a,b)=>b[0]-a[0]).slice(0,topK).forEach(p=>chosen.push(p[1]));
  const firstLo=spec.pairs.findIndex(m=>m.w>lo), firstHi=spec.pairs.findIndex(m=>m.w>hi);
  if(firstLo>=0) chosen.push(firstLo); if(firstHi>=0) chosen.push(firstHi);
  let upper=[]; spec.pairs.forEach((m,i)=>{if(Math.abs(m.w-hi)<0.25) upper.push([spec.ipr[i],i])});
  upper.sort((a,b)=>b[0]-a[0]).slice(0,3).forEach(p=>chosen.push(p[1]));
  spec.pairs.forEach((m,i)=>{if(m.w<lo) chosen.push(i)});
  return Array.from(new Set(chosen)).sort((a,b)=>a-b);
}

function plotLine(id,traces,title){Plotly.newPlot(id,traces,{title,margin:{t:36,l:45,r:15,b:35},paper_bgcolor:'#020617',plot_bgcolor:'#020617',font:{color:'#e5e7eb'}})}

function runSolver(){
  setStatus('Computing background and diagonalizing Hessian...');
  setTimeout(()=>{
    try{
      const gin=+qs('g').value, k1=+qs('k1').value, k2=+qs('k2').value, L=+qs('L').value, Nh=+qs('Nh').value, num=+qs('numModes').value, topK=+qs('topIpr').value;
      const g=qs('sqrt2g').checked?Math.sqrt(2)*gin:gin;
      const bg=mirror(solveBackground(g,k1,k2,L,Nh));
      const sp=Math.sqrt(Math.PI), pc=bg.p1.map((v,i)=>(v+bg.p2[i])/Math.sqrt(2)), ps=bg.p1.map((v,i)=>(v-bg.p2[i])/Math.sqrt(2));
      const vac=eig2(0.5*g*g+2*Math.PI*k1,0.5*g*g,0.5*g*g+2*Math.PI*k2), lo=Math.sqrt(vac[0]), hi=Math.sqrt(vac[1]);
      qs('lowerThr').textContent=lo.toFixed(4); qs('upperThr').textContent=hi.toFixed(4); qs('residual').textContent=bg.res.toExponential(2);
      const spec=spectrum(bg,g,k1,k2,num), sel=chooseModes(spec,lo,hi,topK); qs('selectedModes').textContent=sel.join(', ');
      plotLine('backgroundPhi',[{x:bg.x,y:bg.p1,name:'phi1'},{x:bg.x,y:bg.p2,name:'phi2'}],'Background in phi1, phi2');
      plotLine('backgroundCS',[{x:bg.x,y:pc,name:'phi_c'},{x:bg.x,y:ps,name:'phi_s'}],'Background in phi_c, phi_s');
      const ml=[], mh=[]; for(let i=0;i<bg.x.length;i++){const e=eig2(0.5*g*g+2*Math.PI*k1*Math.cos(2*sp*bg.p1[i]),0.5*g*g,0.5*g*g+2*Math.PI*k2*Math.cos(2*sp*bg.p2[i])); ml.push(e[0]); mh.push(e[1])}
      plotLine('massPlot',[{x:bg.x,y:ml,name:'lower local M2'},{x:bg.x,y:mh,name:'upper local M2'},{x:bg.x,y:bg.x.map(()=>vac[0]),name:'vac lower',line:{dash:'dot'}},{x:bg.x,y:bg.x.map(()=>vac[1]),name:'vac upper',line:{dash:'dot'}}],'Local mass-matrix eigenvalues');
      Plotly.newPlot('spectrumPlot',[{x:spec.pairs.map((_,i)=>i),y:spec.pairs.map(m=>m.w),mode:'markers',marker:{size:spec.ipr.map(v=>7+260*v),color:spec.ipr,colorscale:'Viridis',showscale:true,colorbar:{title:'IPR'}}},{x:[0,spec.pairs.length-1],y:[lo,lo],mode:'lines',name:'lower threshold',line:{dash:'dot'}},{x:[0,spec.pairs.length-1],y:[hi,hi],mode:'lines',name:'upper threshold',line:{dash:'dot'}}],{title:'Spectrum: marker size/color = IPR',margin:{t:36,l:45,r:15,b:35},paper_bgcolor:'#020617',plot_bgcolor:'#020617',font:{color:'#e5e7eb'}});
      const n=bg.x.length-2, xi=bg.x.slice(1,-1), modeTr=[], denTr=[];
      for(const j of sel){const m=spec.pairs[j], arr=[]; let mx=1e-12; for(let r=0;r<n;r++){const val=spec.V[r][m.i]+spec.V[r+n][m.i]; arr.push(val); mx=Math.max(mx,Math.abs(val))} modeTr.push({x:xi,y:arr.map(v=>v/mx),name:'n='+j+', w='+m.w.toFixed(4)}); denTr.push({x:xi,y:spec.dens[j].map(v=>v/Math.max(...spec.dens[j])),name:'n='+j+', IPR='+spec.ipr[j].toFixed(3)});}
      plotLine('modePlot',modeTr,'Selected mode functions, all centered at 0');
      plotLine('densityPlot',denTr,'Selected mode densities');
      let html='<table><tr><th>n</th><th>omega</th><th>omega^2</th><th>IPR</th><th>W(|x|<5)</th><th>class</th></tr>'; spec.pairs.forEach((m,i)=>{if(sel.includes(i)||i<8){const cl=m.w<lo?'below lower':(m.w<hi?'between':'above upper'); html+=`<tr><td>${i}</td><td>${m.w.toFixed(5)}</td><td>${m.w2.toFixed(5)}</td><td>${spec.ipr[i].toFixed(4)}</td><td>${spec.w5[i].toFixed(3)}</td><td>${cl}</td></tr>`}}); html+='</table>'; qs('candidateTable').innerHTML=html;
      setStatus('Done. If a mode looks suspicious, verify it by changing L and checking whether IPR stays O(1).');
    }catch(e){console.error(e); setStatus('Error: '+e.message)}
  },50);
}

qs('runBtn').onclick=runSolver;
