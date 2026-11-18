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
    try {
        return typeof window === 'undefined' ? null : window.localStorage?.getItem(key) ?? null;
    } catch (_error) {
        return null;
    }
}

function _write_local_storage(key, value) {
    try {
        if (typeof window === 'undefined' || !window.localStorage) {
            return false;
        }
        window.localStorage.setItem(key, value);
        return true;
    } catch (_error) {
        return false;
    }
}

function _is_record(value) {
    return value !== null && typeof value === 'object' && !Array.isArray(value);
}

function _is_history_entry(value) {
    return _is_record(value) && typeof value.job_id === 'string' && value.job_id.trim() !== '';
}

function _is_template(value) {
    return _is_record(value) && typeof value.id === 'string' && value.id.trim() !== ''
        && typeof value.name === 'string' && _is_record(value.payload);
}

export function load_task_history() {
    const raw = _read_local_storage(TASK_HISTORY_STORAGE_KEY);
    if (!raw) {
        return [];
    }
    const parsed = _safe_parse(raw, []);
    return Array.isArray(parsed) ? parsed.filter(_is_history_entry).slice(0, TASK_HISTORY_MAX_ITEMS) : [];
}

export function save_task_history(history) {
    const normalized = Array.isArray(history) ? history.filter(_is_history_entry).slice(0, TASK_HISTORY_MAX_ITEMS) : [];
    _write_local_storage(TASK_HISTORY_STORAGE_KEY, JSON.stringify(normalized));
    return normalized;
}

export function upsert_task_history_entry(entry) {
    if (!_is_history_entry(entry)) {
        return load_task_history();
    }
    const history = load_task_history().filter((item) => item.job_id !== entry.job_id);
    history.unshift({ ...entry, updated_at: new Date().toISOString() });
    return save_task_history(history);
}

export function load_task_templates() {
    const raw = _read_local_storage(TASK_TEMPLATE_STORAGE_KEY);
    if (!raw) {
        return get_default_templates();
    }
    const parsed = _safe_parse(raw, []);
    const templates = Array.isArray(parsed) ? parsed.filter(_is_template) : [];
    return templates.length > 0 ? templates : get_default_templates();
}

export function save_task_templates(templates) {
    const normalized = Array.isArray(templates) ? templates.filter(_is_template) : get_default_templates();
    _write_local_storage(TASK_TEMPLATE_STORAGE_KEY, JSON.stringify(normalized));
    return normalized;
}

export function get_default_templates() {
    return JSON.parse(JSON.stringify(DEFAULT_TEMPLATES));
}
