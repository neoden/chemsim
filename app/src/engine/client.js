// Main-thread side of the engine worker.
//   compute(req)  resolves with the result, or with null when a newer request has made this
//                 one obsolete;
//   md(cmd, data) sends a dynamics command; frames come back through onFrame;
//   relax(cmd, data) starts or stops a geometry optimization; frames come back through onRelax.
export function createEngineClient() {
  const worker = new Worker(new URL('./worker.js', import.meta.url), { type: 'module' });
  const waiting = new Map();
  let nextId = 1;
  const settleUpTo = (id, value) => {
    for (const [k, resolve] of waiting) {
      if (k < id) { waiting.delete(k); resolve(null); }
      else if (k === id) { waiting.delete(k); resolve(value); }
    }
  };
  const client = {
    onFrame: null,
    onRelax: null,
    compute(req) {
      const id = nextId++;
      return new Promise((resolve) => {
        waiting.set(id, resolve);
        worker.postMessage({ type: 'compute', ...req, id });
      });
    },
    md(cmd, data = {}) { worker.postMessage({ type: 'md', cmd, ...data }); },
    relax(cmd, data = {}) { worker.postMessage({ type: 'relax', cmd, ...data }); },
  };
  worker.onmessage = ({ data }) => {
    if (data.type === 'md-frame') client.onFrame?.(data);
    else if (data.type === 'relax-frame') client.onRelax?.(data);
    else settleUpTo(data.id, data);
  };
  worker.onerror = (e) => {
    const err = { ok: false, converged: true, error: `The calculation worker crashed: ${e.message || 'unknown error'}` };
    for (const [k, resolve] of waiting) { waiting.delete(k); resolve(err); }
    client.onFrame?.({ type: 'md-frame', running: false, error: err.error, converged: true, fatal: true });
    client.onRelax?.({ type: 'relax-frame', done: true, converged: false, error: err.error, fatal: true });
  };
  return client;
}
