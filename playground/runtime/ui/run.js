const root = document.querySelector('[data-run]');
const status = document.getElementById('runtime-status');
const reason = document.getElementById('runtime-reason');
const restart = document.getElementById('runtime-restart');
const stop = document.getElementById('runtime-stop');
const open = document.getElementById('runtime-open');
const id = root.dataset.run;
const endpoint = '/runtime/api/runs/' + id;
let done = false;
let generation = 0;
let restartRequest = null;
const deadline = performance.now() + Number(root.dataset.deadline) * 1000;
const request = async (url, options = {}) => {
  const response = await fetch(url, {cache:'no-store', credentials:'same-origin',
    signal:AbortSignal.timeout(10000), ...options});
  if (!response.ok) throw Error(response.status === 401 ? 'Authentication is required.' :
    response.status === 403 ? 'This request is not permitted.' :
    response.status === 429 ? 'Application capacity is unavailable.' : 'Runtime is unavailable.');
  return response.json();
};
const mutate = (url, method, body) => request(url, {method,
  headers:{'Content-Type':'application/json', 'X-LCM-Request':'1'},
  ...(body === undefined ? {} : {body:JSON.stringify(body)})});
function render(run) {
  status.textContent = run.state;
  reason.textContent = run.reason || (run.state === 'starting' ? 'Preparing the isolated UI host…' : '');
  const terminal = ['stopped','failed'].includes(run.state);
  restart.hidden = !terminal;
  stop.hidden = terminal;
  open.hidden = run.state !== 'running';
  if (run.state === 'running') {
    const entry = root.dataset.entrypoint;
    const target = root.dataset.entrySurface === 'published' ? new URL(entry, location.origin) :
      new URL('/applications/runs/' + id + entry, location.origin);
    if (root.dataset.entrySurface === 'published') target.searchParams.set('lcm-run', id);
    if (target.origin !== location.origin) throw Error('Invalid registered entry point.');
    open.href = target.href;
    if (root.dataset.automatic === 'true') location.replace(target.href);
  }
  done = terminal || run.state === 'running';
}
restart.addEventListener('click', async () => {
  restart.disabled = true; generation++;
  try {
    const run = await mutate('/runtime/api/runs', 'POST', {
      application:root.dataset.application, request_id:restartRequest ||= crypto.randomUUID()});
    location.assign('/runtime/runs/' + run.id);
  } catch (error) { reason.textContent = error.message; restart.disabled = false; }
});
stop.addEventListener('click', async () => {
  stop.disabled = true; generation++;
  try { render(await mutate(endpoint, 'DELETE')); }
  catch (error) { reason.textContent = error.message; stop.disabled = false; }
});
async function poll() {
  const epoch = generation;
  try {
    const run = await request(endpoint);
    if (epoch !== generation) return;
    render(run);
  } catch (error) {
    if (epoch !== generation) return;
    status.textContent = 'Unavailable'; reason.textContent = error.message;
  }
  if (!done && performance.now() < deadline) setTimeout(poll, 750);
  else if (!done) reason.textContent = 'Status deadline reached. Reload to check this run; no replacement was started.';
}
poll();
