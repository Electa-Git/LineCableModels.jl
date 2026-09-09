// Operator-side fault driver, never browser code or a runtime API. Keeping SSH
// and its process lifecycle here avoids compiling test machinery in a live
// coordinator while its unchanged acknowledgement clock is running.
import {lstat,readFile,writeFile,access} from 'node:fs/promises';
import {isAbsolute,join} from 'node:path';
import {spawn} from 'node:child_process';
import {setTimeout as delay} from 'node:timers/promises';

const [directory,commandFile]=process.argv.slice(2);
const info=await lstat(commandFile);
if(!info.isFile() || info.isSymbolicLink() || (info.mode&0o077) || info.size>8192)
  throw Error('Private operator command file required');
const args=JSON.parse(await readFile(commandFile,'utf8'));
if(!Array.isArray(args) || !args.length || !args.every(v=>typeof v==='string' && !v.includes('\0')) || !isAbsolute(args[0]))
  throw Error('Explicit operator argv required');
let child,closing=false;
process.once('SIGTERM',()=>{closing=true;child?.kill('SIGTERM');});
const exists=path=>access(path).then(()=>true,()=>false);
const exited=(process,seconds)=>new Promise((resolve,reject)=>{
  if(process.exitCode!==null || process.signalCode!==null){resolve({code:process.exitCode,signal:process.signalCode});return;}
  const done=(error,result)=>{clearTimeout(timer);process.off('exit',onExit);process.off('error',onError);error?reject(error):resolve(result);};
  const onExit=(code,signal)=>done(null,{code,signal}),onError=error=>done(error);
  const timer=setTimeout(()=>done(null,null),seconds*1000);
  process.once('exit',onExit);process.once('error',onError);
});
try {
  const deadline=Date.now()+2600_000;
  while(!closing && !await exists(join(directory,'stop-power-worker'))){
    if(Date.now()>=deadline)throw Error('No owned fault request before the fixture deadline');
    await delay(100);
  }
  if(!closing){
    child=spawn(args[0],args.slice(1),{stdio:['ignore','ignore','inherit']});
    const outcome=await exited(child,120);
    if(!outcome || outcome.code!==0)throw Error('Exact remote power-worker stop failed or exceeded its bound');
    await writeFile(join(directory,'power-worker-stopped'),'stopped\n',{mode:0o600});
  }
}finally{
  if(child && child.exitCode===null && child.signalCode===null){
    child.kill('SIGTERM');
    await exited(child,5);
    if(child.exitCode===null && child.signalCode===null)child.kill('SIGKILL');
  }
}
