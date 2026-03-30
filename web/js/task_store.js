const TASK_HISTORY_STORAGE_KEY = 'JRT_TASK_HISTORY_V1';
const TASK_TEMPLATE_STORAGE_KEY = 'JRT_TASK_TEMPLATES_V1';
const TASK_HISTORY_MAX_ITEMS = 50;

const DEFAULT_TEMPLATES = [
    {
        id: 'tmu-smoke',
        name: 'tmu-smoke',
        payload: {
            kind: 'tmu',
            params: {
                T_values: [150.0],
                mu_values: [0.0],
                xi: 0.0,
                max_retries: 0,
            },
        },
    },
    {
        id: 'trho-smoke',
        name: 'trho-smoke',
        payload: {
            kind: 'trho',
            params: {
                T_values: [150.0],
                rho_values: [0.001, 0.002],
                xi: 0.0,
                max_retries: 0,
            },
        },
    },
];

function _safe_parse(raw, fallback) {
    try {
        return JSON.parse(raw);
    } catch (_error) {
        return fallback;
    }
}

function _read_local_storage(key) {
    if (typeof window === 'undefined' || !window.localStorage) {
        return null;
    }
    try {
        return window.localStorage.getItem(key);
    } catch (_error) {
        return null;
    }
}

function _write_local_storage(key, value) {
    if (typeof window === 'undefined' || !window.localStorage) {
        return false;
    }
    try {
        window.localStorage.setItem(key, value);
        return true;
    } catch (_error) {
        return false;
    }
}

export function load_task_history() {
    const raw = _read_local_storage(TASK_HISTORY_STORAGE_KEY);
    if (!raw) {
        return [];
    }
    const parsed = _safe_parse(raw, []);
    return Array.isArray(parsed) ? parsed : [];
}

export function save_task_history(history) {
    const normalized = Array.isArray(history) ? history.slice(0, TASK_HISTORY_MAX_ITEMS) : [];
    _write_local_storage(TASK_HISTORY_STORAGE_KEY, JSON.stringify(normalized));
    return normalized;
}

export function upsert_task_history_entry(entry) {
    if (!entry || !entry.job_id) {
        return load_task_history();
    }
    const history = load_task_history().filter((item) => item.job_id !== entry.job_id);
    history.unshift({ ...entry, updated_at: new Date().toISOString() });
    return save_task_history(history);
}

export function load_task_templates() {
    const raw = _read_local_storage(TASK_TEMPLATE_STORAGE_KEY);
    if (!raw) {
        return DEFAULT_TEMPLATES;
    }
    const parsed = _safe_parse(raw, DEFAULT_TEMPLATES);
    if (!Array.isArray(parsed) || parsed.length === 0) {
        return DEFAULT_TEMPLATES;
    }
    return parsed;
}

export function save_task_templates(templates) {
    const normalized = Array.isArray(templates) ? templates : DEFAULT_TEMPLATES;
    _write_local_storage(TASK_TEMPLATE_STORAGE_KEY, JSON.stringify(normalized));
    return normalized;
}

export function get_default_templates() {
    return DEFAULT_TEMPLATES;
}
