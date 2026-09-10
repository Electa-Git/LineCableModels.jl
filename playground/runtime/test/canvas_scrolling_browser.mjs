// Real CableStudy view: exercise the shared sizing contract through the same
// pointer gestures that used to introduce a scrollbar beside the drawing.
export async function assertCanvasScrolling({read, check, viewport, command, shot, mode}) {
  const settle = () => read('new Promise(resolve => requestAnimationFrame(() => requestAnimationFrame(resolve)))');
  const metrics = () => read(`(() => {
    const view=document.querySelector('.lc-wb-view.is-active');
    const frame=view.querySelector('.lc-viewport-frame');
    const rect=node=>{const r=node.getBoundingClientRect();return {top:r.top,bottom:r.bottom,left:r.left,right:r.right,width:r.width,height:r.height};};
    const body=frame.querySelector('.lc-viewport-body');
    const circle=view.querySelector('.lc-study-sheath');
    const inners=[...view.querySelectorAll('.lc-wb-split-mount,.lc-wb-split,.lc-wb-split-region,.lc-viewport-body,.lc-viewport-content')];
    return {view:rect(view),frame:rect(frame),body:rect(body),circle:rect(circle),
      scroll:view.scrollHeight-view.clientHeight, horizontal:view.scrollWidth-view.clientWidth,
      innerScrollers:inners.filter(n=>['auto','scroll'].includes(getComputedStyle(n).overflowY) && n.scrollHeight>n.clientHeight+1).map(n=>n.className),
      clipping:inners.filter(n=>['hidden','clip'].includes(getComputedStyle(n).overflowY) && n.scrollHeight>n.clientHeight+1).map(n=>n.className)};
  })()`);
  const fitted = async label => {
    const m=await metrics();
    await check(`${m.innerScrollers.length === 0 && m.clipping.length === 0 && m.horizontal <= 1}`,
      label+' introduced an inner scroll/clipped surface: '+JSON.stringify(m));
    await check(`${Math.abs(m.circle.width-m.circle.height)<1 && m.circle.width>30 &&
      m.circle.top>=m.body.top && m.circle.bottom<=m.body.bottom && m.circle.left>=m.body.left && m.circle.right<=m.body.right}`,
      label+' distorted or clipped the drawing: '+JSON.stringify(m));
    return m;
  };
  const dock = async state => {
    await read(`(() => {if(document.querySelector('.lc-wb-shell').dataset.dockState!==${JSON.stringify(state)})
      document.querySelector('.lc-wb-dock-toggle').click();})()`);
    await settle();
  };
  await viewport(1920,1047); await dock('collapsed'); await settle();
  await read(`window.keptGeometry=document.querySelector('.lc-wb-view.is-active .lc-study-plot')`);
  const original=await fitted(mode+' initial construction');
  await check(`${original.scroll<=1 && original.frame.bottom<=original.view.bottom}`,
    mode+' fitting construction must not scroll: '+JSON.stringify(original));
  for (const delta of [-240,480,-240]) {
    const grip=await read(`(() => {const r=document.querySelector('.lc-wb-view.is-active .lc-wb-splitter-handle').getBoundingClientRect();return {x:r.x+r.width/2,y:r.y+r.height/2};})()`);
    await command('Input.dispatchMouseEvent',{type:'mouseMoved',...grip});
    await command('Input.dispatchMouseEvent',{type:'mousePressed',...grip,button:'left',buttons:1,clickCount:1});
    for (let step=1;step<=8;step++) {
      await command('Input.dispatchMouseEvent',{type:'mouseMoved',x:grip.x+delta*step/8,y:grip.y,button:'left',buttons:1});
      await settle();
      const m=await fitted(mode+' construction drag '+delta+' step '+step);
      await check(`${Math.abs(m.frame.height-original.frame.height)<=1 && m.scroll<=1}`,
        mode+' widening the drawing changed height or created overflow: '+JSON.stringify(m));
    }
    await command('Input.dispatchMouseEvent',{type:'mouseReleased',x:grip.x+delta,y:grip.y,button:'left',clickCount:1});
  }
  await check(`keptGeometry===document.querySelector('.lc-wb-view.is-active .lc-study-plot')`, mode+' split drag remounted geometry');
  // The optional frame cap is independent of pane width and includes the header.
  await read(`keptGeometry.closest('.lc-viewport-frame').style.maxHeight='320px'`); await settle();
  const capped=await fitted(mode+' capped construction');
  await check(`${Math.abs(capped.frame.height-320)<=1}`, mode+' frame maximum height ignored');
  await read(`keptGeometry.closest('.lc-viewport-frame').style.removeProperty('max-height')`);
  await dock('expanded'); await settle();
  const expanded=await fitted(mode+' expanded diagnostics');
  await check(`${expanded.frame.height<original.frame.height && expanded.scroll<=1}`, mode+' dock expansion did not refit the frame');
  await shot('construction-'+mode);
  await dock('collapsed');
  for (const [width,height] of [[1920,360],[390,800]]) {
    await viewport(width,height); await settle();
    const m=await fitted(mode+' construction '+width+'x'+height);
    await check(`${m.scroll>1}`, mode+' overflowing content did not reach the outer canvas');
    await check(`(() => {const view=document.querySelector('.lc-wb-view.is-active'),workspace=document.querySelector('.lc-wb-workspace');
      return getComputedStyle(view).overflowY==='auto' && Math.abs(view.getBoundingClientRect().right-workspace.getBoundingClientRect().right)<1;
    })()`, mode+' scrolling is not at the outer canvas edge');
    await read(`document.querySelector('.lc-wb-view.is-active .lc-property-grid').scrollIntoView({block:'end'})`); await settle();
    await check(`(() => {const view=document.querySelector('.lc-wb-view.is-active').getBoundingClientRect(),
      last=document.querySelector('.lc-wb-view.is-active .lc-property-grid').getBoundingClientRect();
      return last.bottom<=view.bottom+1 && last.top>=view.top && document.documentElement.scrollWidth<=innerWidth+1;
    })()`, mode+' outer scrolling cannot reach the complete inputs at '+width);
    await shot('construction-'+mode+'-'+width);
  }
  await viewport(1440,900); await dock('expanded');
  await read(`document.querySelector('.lc-wb-view.is-active').scrollTop=0`); await settle();
}
