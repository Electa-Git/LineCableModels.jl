// Published widget embeds share one run prefix. This script never allocates a
// process: the catalogue launcher owns that explicit action.
const frames=[...document.querySelectorAll('iframe[data-lc-published-src]')];
if(frames.length) {
  const query=new URLSearchParams(location.search);
  const run=query.get('lcm-run');
  const validRun=run && /^[a-f0-9]{8}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{4}-[a-f0-9]{12}$/.test(run);
  let runtime=false;
  try {
    const response=await fetch('/runtime/api/capabilities',{cache:'no-store',signal:AbortSignal.timeout(3000)});
    runtime=response.ok && (await response.json()).ui_hosts===true;
  } catch (_) { /* The developer publisher retains its standalone routes. */ }
  for(const frame of frames) {
    const path=frame.dataset.lcPublishedSrc;
    const safe=path?.startsWith('/') && !path.startsWith('//') && !/[\\]|(?:%2f|%5c|%2e)|\/\.\.?(?:\/|$)/i.test(path);
    // Legacy v1 jobs never enter an authenticated v2 run by fallback.
    const unsupported=runtime && path==='/widgets/job-panel';
    if(!safe || (query.has('lcm-run') && !validRun) || (runtime && !validRun) || unsupported) {
      frame.hidden=true;
      const notice=document.createElement('p'); notice.className='lc-catalogue-status';
      notice.textContent=unsupported ? 'The legacy unassigned job panel is not exposed in an owned runtime. Assigned controls are installed separately.' :
        'Interactive preview requires an explicitly started toolkit run.';
      frame.after(notice);
    } else {
      frame.src=validRun ? '/applications/runs/'+run+path : path;
    }
  }
  if(runtime && !validRun && !query.has('lcm-run')) {
    const link=document.createElement('a'); link.href='/dev/';link.textContent='Start the interactive toolkit gallery →';
    document.querySelector('main')?.prepend(link);
  }
}
