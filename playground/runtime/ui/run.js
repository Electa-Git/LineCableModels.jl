const root = document.querySelector('[data-run]');
const status = document.getElementById('runtime-status');
const reason = document.getElementById('runtime-reason');
const restart = document.getElementById('runtime-restart');
const stop = document.getElementById('runtime-stop');
const open = document.getElementById('runtime-open');
// A gateway may still serve its previously loaded Julia markup while adopting
// new static assets. Enhance that markup without requiring a service restart.
const elapsed = document.getElementById('runtime-elapsed') || document.createElement('p');
if (!elapsed.isConnected) {
  elapsed.id = 'runtime-elapsed';
  elapsed.className = 'lc-runtime-hint';
  elapsed.setAttribute('aria-live', 'off');
  reason.after(elapsed);
}
status.classList.add('lc-activity-status');
const id = root.dataset.run;
const endpoint = '/runtime/api/runs/' + id;
let done = false;
let generation = 0;
let restartRequest = null;
let pollTimer;
let elapsedTimer;
const openedAt = performance.now();
const deadline = openedAt + Number(root.dataset.deadline) * 1000;
function activity(label, busy, tone = '', waitingForHost = busy) {
  status.textContent = label;
  status.dataset.busy = String(busy);
  status.dataset.tone = tone;
  clearInterval(elapsedTimer);
  elapsed.hidden = !waitingForHost;
  if (waitingForHost) {
    const tick = () => {
      const seconds = Math.floor((performance.now() - openedAt) / 1000);
      elapsed.textContent = `Waiting for application readiness · ${seconds} s on this page. First startup may take longer while Julia loads and compiles.`;
    };
    tick(); elapsedTimer = setInterval(tick, 1000);
  }
}
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
  const busy = ['reserved', 'starting'].includes(run.state);
  activity(run.state, busy, run.state === 'failed' ? 'danger' : '');
  reason.textContent = run.reason || (busy ? 'Preparing the isolated UI host…' : '');
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
  clearTimeout(pollTimer);
  activity('Requesting a clean run…', true, '', false);
  try {
    const run = await mutate('/runtime/api/runs', 'POST', {
      application:root.dataset.application, request_id:restartRequest ||= crypto.randomUUID()});
    location.assign('/runtime/runs/' + run.id);
  } catch (error) {
    activity('Restart not confirmed', false, 'warning');
    reason.textContent = error.message; restart.disabled = false;
  }
});
stop.addEventListener('click', async () => {
  stop.disabled = true; generation++;
  clearTimeout(pollTimer);
  activity('Stopping run…', true, '', false);
  try { render(await mutate(endpoint, 'DELETE')); }
  catch (error) {
    activity('Stop not confirmed', false, 'warning');
    reason.textContent = error.message; stop.disabled = false;
    schedulePoll();
  }
});
function schedulePoll() {
  if (done) return;
  if (performance.now() < deadline) pollTimer = setTimeout(poll, 750);
  else {
    activity('Status check paused', false, 'warning');
    reason.textContent = 'Status deadline reached. Reload to check this run; no replacement was started.';
  }
}
async function poll() {
  const epoch = generation;
  try {
    const run = await request(endpoint);
    if (epoch !== generation) return;
    render(run);
  } catch (error) {
    if (epoch !== generation) return;
    activity('Status unavailable', false, 'warning');
    reason.textContent = error.message + ' Retrying status checks; no replacement run is being started.';
  }
  schedulePoll();
}
activity(status.textContent, ['reserved', 'starting'].includes(status.textContent.trim()));
addEventListener('pagehide', () => {
  generation++; clearTimeout(pollTimer); clearInterval(elapsedTimer);
});
addEventListener('pageshow', event => { if (event.persisted) poll(); });
poll();
