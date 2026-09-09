// One passive catalogue is emitted from ApplicationCatalogue.entries() by lcm.
// Reading a page or changing the selection never starts a run or a worker.
async function mountCatalogue(root) {
  const response = await fetch('/assets/application-catalogue.json', {cache:'no-cache'});
  if (!response.ok) throw Error('Catalogue is unavailable. Rebuild with lcm playground build.');
  const records = await response.json();
  const entries = records.filter(entry => entry.kind === root.dataset.kind &&
    (root.dataset.visibility === 'all' || entry.visibility === (root.dataset.visibility || 'public')));
  const hint = document.createElement('p'); hint.className='lc-catalogue-status'; hint.role='status';
  if (!entries.length) {
    hint.textContent='No public application is registered in this category yet. The developer foundation remains available below.';
    root.append(hint); return;
  }
  const chooser=document.createElement('div'); chooser.className='lc-catalogue-chooser';
  const label=document.createElement('label'); label.append(document.createTextNode('Choose ' + root.dataset.kind));
  const select=document.createElement('select'); select.className='lc-control-select lc-form-control';
  for (const entry of entries) {
    const option=document.createElement('option'); option.value=entry.id; option.textContent=entry.title;
    select.append(option);
  }
  label.append(select); chooser.append(label);
  const actions=document.createElement('div'); actions.className='lc-catalogue-actions';
  const read=document.createElement('a'); read.textContent='Open static slides';
  const launch=document.createElement('button'); launch.className='lc-button lc-button-secondary'; launch.type='button'; launch.textContent='Start isolated run'; launch.disabled=true;
  actions.append(read,launch); chooser.append(actions);
  const description=document.createElement('p');
  const roles=document.createElement('p'); roles.className='lc-catalogue-roles';
  root.append(chooser,description,roles,hint);
  let available=false, installed=new Set(), requestId=null, busy=false;
  const selected=()=>entries.find(entry=>entry.id===select.value);
  function update() {
    const entry=selected();
    description.textContent=entry.description;
    roles.textContent=entry.requirements.length ? 'Declared roles: ' + entry.requirements.map(r=>r.role).join(' · ') + '. Preparation is explicit.' :
      'Local UI demonstration. No scientific worker is required.';
    read.hidden=entry.kind!=='presentation'; read.href=entry.entrypoint;
    launch.disabled=!available || !installed.has(entry.id) || busy;
    hint.textContent=available && !installed.has(entry.id) ? 'Static content available · this application’s live implementation is not installed.' :
      available ? 'UI host available · scientific resources are not prepared by opening this page.' :
      'Static publishing mode · runtime launcher is unavailable. Slides and developer documentation remain usable.';
  }
  select.addEventListener('change',()=>{requestId=null; update();});
  launch.addEventListener('click',async()=>{
    if (busy) return;
    busy=true; select.disabled=true; launch.disabled=true;
    requestId ||= crypto.randomUUID();
    try {
      const response=await fetch('/runtime/api/runs', {method:'POST', credentials:'same-origin',
        headers:{'Content-Type':'application/json','X-LCM-Request':'1'},
        body:JSON.stringify({application:selected().id,request_id:requestId}), signal:AbortSignal.timeout(15000)});
      if (!response.ok) throw Error(response.status===401 ? 'Authenticate to start a private run.' :
        response.status===429 ? 'Application capacity is unavailable.' : 'This run could not be started.');
      const run=await response.json();
      if (!/^[a-f0-9-]{36}$/.test(run.id)) throw Error('Invalid runtime response.');
      location.assign('/runtime/runs/' + run.id);
    } catch(error) {hint.textContent=error.message; busy=false; select.disabled=false; launch.disabled=!available;}
  });
  update();
  try {
    const capability=await fetch('/runtime/api/capabilities', {signal:AbortSignal.timeout(3000),cache:'no-store'});
    if (capability.ok) {
      const data=await capability.json();
      available=data.ui_hosts===true; installed=new Set(data.applications || []);
    }
  } catch (_) {available=false;}
  update();
}
for (const root of document.querySelectorAll('[data-lcm-application-catalogue]')) {
  mountCatalogue(root).catch(error=>{const notice=document.createElement('p');notice.textContent=error.message;root.append(notice);});
}
