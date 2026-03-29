const DEFAULT_API_BASE_URL = 'http://localhost:8080';

function normalize_base_url(raw) {
    const value = typeof raw === 'string' ? raw.trim() : '';
    if (!value) {
        return DEFAULT_API_BASE_URL;
    }
    return value.replace(/\/+$/, '');
}

export function get_api_base_url() {
    if (typeof globalThis !== 'undefined') {
        const injected = globalThis.__JRT_API_BASE_URL__;
        if (typeof injected === 'string' && injected.trim()) {
            return normalize_base_url(injected);
        }
    }

    if (typeof window !== 'undefined') {
        const urlParam = new URLSearchParams(window.location.search).get('api_base_url');
        if (urlParam && urlParam.trim()) {
            return normalize_base_url(urlParam);
        }

        try {
            const storageValue = window.localStorage.getItem('JRT_API_BASE_URL');
            if (storageValue && storageValue.trim()) {
                return normalize_base_url(storageValue);
            }
        } catch (_error) {
            return DEFAULT_API_BASE_URL;
        }
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

export { DEFAULT_API_BASE_URL };
