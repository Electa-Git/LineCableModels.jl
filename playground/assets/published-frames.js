// Published widget embeds share one run prefix. This script never allocates a
// process: the catalogue launcher owns that explicit action.
const frames=[...document.querySelectorAll('iframe[data-lc-published-src]')];
if(frames.length) {
  const styleURL=new URL('./published-frames.css', import.meta.url).href;
  if (![...document.querySelectorAll('link[rel="stylesheet"]')].some(link=>link.href===styleURL)) {
    const style=document.createElement('link'); style.rel='stylesheet'; style.href=styleURL;
    document.head.append(style);
  }
  // Older open documents can adopt the new shared embed without a server
  // restart. New shortcode output already has this exact structure.
  for (const frame of frames) {
    if (frame.parentElement.classList.contains('lc-published-viewport')) continue;
    const viewport=document.createElement('div');
    viewport.className='lc-published-viewport'; viewport.dataset.lcPreviewState='pending';
    viewport.style.setProperty('--lc-widget-height', frame.style.getPropertyValue('--lc-widget-height') || '20rem');
    const placeholder=document.createElement('div'); placeholder.className='lc-published-placeholder';
    placeholder.setAttribute('role','status');
    const title=document.createElement('strong'); title.textContent=frame.title;
    const message=document.createElement('p'); message.textContent='Interactive preview is not connected.';
    placeholder.append(title,message); frame.before(viewport); viewport.append(frame,placeholder);
  }
  const query=new URLSearchParams(location.search);
  const run=query.get('lcm-run');
  const validRun=run && /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/.test(run);
  let runtime=null;
  try {
    const response=await fetch('/runtime/api/capabilities',{cache:'no-store',signal:AbortSignal.timeout(3000)});
    // A standalone publisher has no capabilities endpoint. A failed gateway
    // check is not permission to fall through to a different route namespace.
    if (response.status===404) runtime=false;
    else if (response.ok) {
      const capabilities=await response.json();
      if (typeof capabilities.ui_hosts==='boolean') runtime=capabilities.ui_hosts;
    }
  } catch (_) { /* Unknown availability stays in the owned inactive surface. */ }
  for(const frame of frames) {
    const viewport=frame.parentElement;
    const placeholder=viewport.querySelector('.lc-published-placeholder');
    const path=frame.dataset.lcPublishedSrc;
    const safe=path?.startsWith('/') && !path.startsWith('//') && !/[\\]|(?:%2f|%5c|%2e)|\/\.\.?(?:\/|$)/i.test(path);
    // Legacy v1 jobs never enter an authenticated v2 run by fallback.
    const unsupported=runtime && path==='/widgets/job-panel';
    if(!safe || runtime===null || (query.has('lcm-run') && !validRun) || (runtime && !validRun) || unsupported) {
      frame.hidden=true;
      placeholder.hidden=false;
      viewport.dataset.lcPreviewState='unavailable';
      placeholder.querySelector('p').textContent=runtime===null ? 'Interactive preview service is unavailable. Reload this page to retry; no application was started.' :
        unsupported ? 'The legacy unassigned job panel is not exposed in an owned runtime. Assigned controls are installed separately.' :
        'Interactive preview requires an explicitly started toolkit run.';
    } else {
      frame.hidden=false;
      placeholder.hidden=true;
      viewport.dataset.lcPreviewState='live';
      // Establish the real rectangle before assigning src. A display:none
      // frame has no distance from the viewport and can bypass lazy loading,
      // eagerly booting every gallery app while blocking the page's load.
      frame.getBoundingClientRect();
      frame.src=validRun ? '/applications/runs/'+run+path : path;
    }
  }
  if(runtime && !validRun && !query.has('lcm-run')) {
    const link=document.createElement('a'); link.href='/dev/';link.textContent='Start the interactive toolkit gallery →';
    // A launcher is document content, not a stray inline node above the h1.
    const notice=document.createElement('p'); notice.className='lc-catalogue-launcher'; notice.append(link);
    const main=document.querySelector('main');
    const title=main?.querySelector(':scope > .quarto-title-block');
    if(title) title.after(notice); else main?.prepend(notice);
  }
}
