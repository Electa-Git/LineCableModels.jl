import assert from 'node:assert/strict';

// Uses the real gateway HTML/assets and an explicit browser-side HTTP fixture.
// No stop/restart request in this test may reach a real application host.
export async function assertRuntimeActivity({command, read, wait, base, run, shot}) {
  const endpoint = '/runtime/api/runs/' + run.id;
  const {identifier} = await command('Page.addScriptToEvaluateOnNewDocument', {source: `
    window.runEvidence = {state:'starting', reason:'Preparing the isolated UI host…'};
    window.runtimeDown = false; window.stopFails = false; window.runtimeRequests = [];
    const realFetch = window.fetch;
    window.fetch = async (url, options = {}) => {
      const path = new URL(url, location.href).pathname;
      if (!path.startsWith('/runtime/api/runs')) return realFetch(url, options);
      runtimeRequests.push(options.method || 'GET');
      if (path !== ${JSON.stringify(endpoint)}) throw Error('Unexpected runtime mutation in activity fixture');
      if (runtimeDown || (options.method === 'DELETE' && stopFails)) return new Response('{}', {status:503});
      if (options.method === 'DELETE') runEvidence = {state:'stopped', reason:'Stopped by owner'};
      return new Response(JSON.stringify(runEvidence), {headers:{'Content-Type':'application/json'}});
    };
    const observer = new MutationObserver(() => {
      const root = document.querySelector('[data-run]');
      if (!root) return;
      root.dataset.automatic = 'false';
      if (location.hash === '#deadline') root.dataset.deadline = '0.1';
      if (location.hash === '#legacy') {
        const elapsed=document.querySelector('#runtime-elapsed');
        if (!elapsed) return;
        elapsed.remove();
        document.querySelector('#runtime-status').className='';
      }
      observer.disconnect();
    });
    observer.observe(document, {childList:true, subtree:true});
  `});
  let navigation = 0;
  const navigate = async (suffix = '') => {
    const url=base + '/runtime/runs/' + run.id + '?activity=' + (++navigation) + suffix;
    await command('Page.navigate', {url});
    await wait(`location.href===${JSON.stringify(url)} && document.readyState==='complete' &&
      (document.querySelector('#runtime-status')?.textContent==='starting' ||
      document.querySelector('#runtime-status')?.textContent==='Status check paused')`, 'Runtime activity did not mount');
  };
  const setState = async state => {
    await read(`runEvidence={state:${JSON.stringify(state)},reason:''}`);
    await wait(`document.querySelector('#runtime-status').textContent===${JSON.stringify(state)}`, 'Missing runtime state ' + state);
  };
  const busy = () => read(`document.querySelector('#runtime-status').dataset.busy`);
  try {
    for (const theme of ['dark','light']) {
      await navigate();
      await read(`LineCableModelsTheme.select('${theme}')`);
      assert.equal(await busy(), 'true');
      assert.equal(await read(`(() => {
        const node=document.querySelector('#runtime-status'), s=getComputedStyle(node,'::before');
        return node.getAttribute('role')==='status' && s.animationName==='lc-activity-spin' &&
          s.borderTopColor!==s.borderBottomColor;
      })()`), true);
      await wait(`!document.querySelector('#runtime-elapsed').textContent.includes('· 0 s')`, 'Elapsed indication did not advance');
      for (const width of [1440,390]) {
        await command('Emulation.setDeviceMetricsOverride', {width,height:900,deviceScaleFactor:1,mobile:false});
        assert.equal(await read('document.documentElement.scrollWidth <= innerWidth'),true);
        await shot('startup-'+theme+'-'+width);
      }
      await command('Emulation.setEmulatedMedia', {features:[{name:'prefers-reduced-motion',value:'reduce'}]});
      assert.equal(await read(`getComputedStyle(document.querySelector('#runtime-status'),'::before').animationName`),'none');
      await command('Emulation.setEmulatedMedia', {features:[]});
      await read('runtimeDown=true');
      await wait(`document.querySelector('#runtime-status').textContent==='Status unavailable'`, 'Status loss not visible');
      assert.equal(await busy(),'false');
      await read('runtimeDown=false');
      await wait(`document.querySelector('#runtime-status').textContent==='starting'`, 'Status did not recover');
      assert.equal(await busy(),'true');
      await setState('running');
      assert.equal(await busy(),'false');
      assert.equal(await read(`document.querySelector('#runtime-elapsed').hidden && !document.querySelector('#runtime-open').hidden`),true);
      // Terminal failure/stop are distinct from an unknown/unavailable status.
      for (const state of ['failed','stopped']) {
        await navigate(); await setState(state);
        assert.equal(await busy(),'false');
        assert.equal(await read(`!document.querySelector('#runtime-restart').hidden && document.querySelector('#runtime-stop').hidden`),true);
      }
    }
    await navigate();
    await read(`stopFails=true; document.querySelector('#runtime-stop').click()`);
    await wait(`document.querySelector('#runtime-status').textContent==='Stop not confirmed'`, 'Stop failure not visible');
    assert.equal(await busy(),'false');
    await wait(`document.querySelector('#runtime-status').textContent==='starting'`, 'Failed stop permanently stalled polling');
    await read(`stopFails=false; document.querySelector('#runtime-stop').click()`);
    await wait(`document.querySelector('#runtime-status').textContent==='stopped'`, 'Successful stop not visible');
    assert.equal(await busy(),'false');
    await navigate('#deadline');
    await wait(`document.querySelector('#runtime-status').textContent==='Status check paused'`, 'Deadline still implies active startup');
    assert.equal(await busy(),'false');
    assert.equal(await read(`runtimeRequests.every(method=>method==='GET')`),true);
    await navigate('#legacy');
    assert.equal(await busy(),'true');
    assert.equal(await read(`document.querySelectorAll('#runtime-elapsed').length===1 &&
      document.querySelector('#runtime-status').classList.contains('lc-activity-status')`),true);
    console.log('Runtime activity: startup, elapsed, readiness, failure, stop/recovery, deadline, reduced motion and narrow layout; both themes');
  } finally {
    await command('Page.removeScriptToEvaluateOnNewDocument',{identifier});
    await command('Emulation.setEmulatedMedia',{features:[]});
  }
}
