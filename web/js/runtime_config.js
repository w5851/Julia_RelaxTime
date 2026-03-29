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

export function get_api_base_url() {
    if (typeof globalThis !== 'undefined') {
        const injectedProfile = globalThis.__JRT_DEPLOY_PROFILE__;
        if (typeof injectedProfile === 'string' && injectedProfile.trim()) {
            return resolve_profile_api_base_url(injectedProfile);
        }

        const injected = globalThis.__JRT_API_BASE_URL__;
        if (typeof injected === 'string' && injected.trim()) {
            return normalize_base_url(injected);
        }
    }

    if (typeof window !== 'undefined') {
        const profileParam = new URLSearchParams(window.location.search).get('deploy_profile');
        if (profileParam && profileParam.trim()) {
            return resolve_profile_api_base_url(profileParam);
        }

        const urlParam = new URLSearchParams(window.location.search).get('api_base_url');
        if (urlParam && urlParam.trim()) {
            return normalize_base_url(urlParam);
        }

        try {
            const profileStorageValue = window.localStorage.getItem('JRT_DEPLOY_PROFILE');
            if (profileStorageValue && profileStorageValue.trim()) {
                return resolve_profile_api_base_url(profileStorageValue);
            }

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
