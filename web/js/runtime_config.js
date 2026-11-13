const DEFAULT_API_BASE_URL = 'http://localhost:8080';
const DEFAULT_DEPLOY_PROFILE = 'localhost';
const PROFILE_API_BASE_URLS = {
    localhost: 'http://localhost:8080',
    staging: 'https://staging.jrt.local',
    remote: 'https://api.jrt.example.com',
};

function normalize_deploy_profile(raw) {
    const value = typeof raw === 'string' ? raw.trim().toLowerCase() : '';
    if (!value) {
        return DEFAULT_DEPLOY_PROFILE;
    }
    if (Object.prototype.hasOwnProperty.call(PROFILE_API_BASE_URLS, value)) {
        return value;
    }
    return DEFAULT_DEPLOY_PROFILE;
}

export function resolve_profile_api_base_url(profile) {
    const normalized = normalize_deploy_profile(profile);
    return PROFILE_API_BASE_URLS[normalized] || DEFAULT_API_BASE_URL;
}

function normalize_base_url(raw) {
    const value = typeof raw === 'string' ? raw.trim() : '';
    if (!value) {
        return DEFAULT_API_BASE_URL;
    }
    return value.replace(/\/+$/, '');
}

function stored_value(key) {
    try {
        return typeof window !== 'undefined' ? window.localStorage.getItem(key) : null;
    } catch (_error) {
        return null;
    }
}

function query_value(key) {
    return typeof window !== 'undefined'
        ? new URLSearchParams(window.location.search).get(key) : null;
}

function nonempty_string(value) {
    return typeof value === 'string' && value.trim() ? value : null;
}

export function get_api_base_url() {
    const manual = nonempty_string(globalThis.__JRT_API_BASE_URL__)
        || nonempty_string(query_value('api_base_url'))
        || nonempty_string(stored_value('JRT_API_BASE_URL'));
    if (manual) {
        return normalize_base_url(manual);
    }
    const profile = nonempty_string(globalThis.__JRT_DEPLOY_PROFILE__)
        || nonempty_string(query_value('deploy_profile'))
        || nonempty_string(stored_value('JRT_DEPLOY_PROFILE'));
    if (profile) {
        return resolve_profile_api_base_url(profile);
    }
    if (typeof window !== 'undefined' && /^https?:$/.test(window.location.protocol)) {
        return window.location.origin;
    }
    return DEFAULT_API_BASE_URL;
}

export function set_api_base_url(value) {
    const normalized = normalize_base_url(value);
    if (typeof window !== 'undefined') {
        try {
            window.localStorage.setItem('JRT_API_BASE_URL', normalized);
        } catch (_error) {
            return normalized;
        }
    }
    return normalized;
}

export function get_deploy_profile() {
    if (typeof globalThis !== 'undefined') {
        const injected = globalThis.__JRT_DEPLOY_PROFILE__;
        if (typeof injected === 'string' && injected.trim()) {
            return normalize_deploy_profile(injected);
        }
    }

    if (typeof window !== 'undefined') {
        const urlParam = new URLSearchParams(window.location.search).get('deploy_profile');
        if (urlParam && urlParam.trim()) {
            return normalize_deploy_profile(urlParam);
        }
        try {
            const storageValue = window.localStorage.getItem('JRT_DEPLOY_PROFILE');
            if (storageValue && storageValue.trim()) {
                return normalize_deploy_profile(storageValue);
            }
        } catch (_error) {
            return DEFAULT_DEPLOY_PROFILE;
        }
    }

    return DEFAULT_DEPLOY_PROFILE;
}

export function set_deploy_profile(value) {
    const normalized = normalize_deploy_profile(value);
    if (typeof window !== 'undefined') {
        try {
            window.localStorage.setItem('JRT_DEPLOY_PROFILE', normalized);
        } catch (_error) {
            return normalized;
        }
    }
    return normalized;
}

export { DEFAULT_API_BASE_URL, DEFAULT_DEPLOY_PROFILE };
