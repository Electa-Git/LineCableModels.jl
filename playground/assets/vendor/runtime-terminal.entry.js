// Pinned, locally served terminal renderer only. The owned component supplies
// authorization, bounded transport, theme projection and lifetime management.
import { Terminal } from '@xterm/xterm';
import { FitAddon } from '@xterm/addon-fit';
import '@xterm/xterm/css/xterm.css';

globalThis.LineCableModelsTerminalVendor ??= Object.freeze({ Terminal, FitAddon });
