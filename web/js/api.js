/**
 * API通信模块
 * 负责与Julia HTTP服务器通信
 */

import { get_api_base_url } from './runtime_config.js';

const DEFAULT_TIMEOUT_MS = 10000;
const DEFAULT_RETRIES = 1;

export function build_api_url(path) {
    const base = get_api_base_url();
    const normalizedPath = path.startsWith('/') ? path : `/${path}`;
    return `${base}${normalizedPath}`;
}

function is_retryable_error(error) {
    return error && (error.name === 'AbortError' || error.name === 'TypeError');
}

async function request_json(path, options = {}) {
    const {
        timeoutMs = DEFAULT_TIMEOUT_MS,
        retries = DEFAULT_RETRIES,
        ...fetchOptions
    } = options;

    let attempt = 0;
    while (true) {
        const controller = new AbortController();
        const timeoutId = setTimeout(() => controller.abort(), timeoutMs);
        try {
            const response = await fetch(build_api_url(path), {
                ...fetchOptions,
                signal: controller.signal,
                headers: {
                    'Content-Type': 'application/json',
                    ...(fetchOptions.headers || {}),
                },
            });
            clearTimeout(timeoutId);

            const payload = await response.json().catch(() => ({}));
            if (!response.ok) {
                throw normalize_api_error({
                    status: response.status,
                    payload,
                    message: payload.message || payload.error || `HTTP error! status: ${response.status}`,
                });
            }
            return payload;
        } catch (error) {
            clearTimeout(timeoutId);
            if (attempt < retries && is_retryable_error(error)) {
                attempt += 1;
                continue;
            }
            if (error && error.isNormalizedApiError) {
                throw error;
            }
            throw normalize_api_error({
                status: null,
                payload: {},
                message: error?.message || '请求失败',
                originalError: error,
            });
        }
    }
}

export function normalize_api_error({ status, payload = {}, message = '', originalError = null }) {
    const code = payload.code || payload.error_code || (status === 429 ? 'QUEUE_FULL' : 'API_ERROR');
    const normalized = new Error(message || payload.message || payload.error || '请求失败');
    normalized.code = code;
    normalized.status = status;
    normalized.diagnostics = payload.diagnostics || null;
    normalized.payload = payload;
    normalized.originalError = originalError;
    normalized.isNormalizedApiError = true;
    return normalized;
}

function _is_number(value) {
    return Number.isFinite(Number(value));
}

export function validate_scan_request(request) {
    const errors = [];
    const kind = String(request?.kind || '').toLowerCase();
    if (!(kind === 'tmu' || kind === 'trho')) {
        errors.push('kind 必须是 tmu 或 trho');
    }

    const params = request?.params || {};
    const hasXi = params.xi !== undefined && params.xi !== null && params.xi !== '';
    const hasXiValues = Array.isArray(params.xi_values) && params.xi_values.length > 0;
    const hasXiGrid = !!params.xi_grid;
    const xiStrategies = [hasXi, hasXiValues, hasXiGrid].filter(Boolean).length;
    if (xiStrategies !== 1) {
        errors.push('xi 策略必须且只能选择一种：xi / xi_values / xi_grid');
    }

    if (hasXi && !_is_number(params.xi)) {
        errors.push('xi 必须是数字');
    }

    if (hasXiValues) {
        const allNumeric = params.xi_values.every(v => _is_number(v));
        if (!allNumeric) {
            errors.push('xi_values 必须全部是数字');
        }
    }

    if (hasXiGrid) {
        const start = params.xi_grid.start;
        const stop = params.xi_grid.stop;
        const step = params.xi_grid.step;
        if (!(_is_number(start) && _is_number(stop) && _is_number(step))) {
            errors.push('xi_grid.start/stop/step 必须是数字');
        } else if (Number(step) <= 0) {
            errors.push('xi_grid.step 必须大于 0');
        } else if (Number(stop) < Number(start)) {
            errors.push('xi_grid.stop 必须大于等于 xi_grid.start');
        }
    }

    if (!Array.isArray(params.T_values) || params.T_values.length === 0 || !params.T_values.every(v => _is_number(v))) {
        errors.push('T_values 必须是非空数字数组');
    }

    if (kind === 'tmu') {
        if (!Array.isArray(params.mu_values) || params.mu_values.length === 0 || !params.mu_values.every(v => _is_number(v))) {
            errors.push('tmu 模式下 mu_values 必须是非空数字数组');
        }
    }

    if (kind === 'trho') {
        if (!Array.isArray(params.rho_values) || params.rho_values.length === 0 || !params.rho_values.every(v => _is_number(v))) {
            errors.push('trho 模式下 rho_values 必须是非空数字数组');
        }
    }

    return {
        valid: errors.length === 0,
        errors,
    };
}

export class API {
    /**
     * 检查服务器健康状态
     * @returns {Promise<boolean>} 服务器是否在线
     */
    static async checkHealth() {
        const url = build_api_url('/health');
        console.log(`[API] Checking health at ${url}`);

        const maxAttempts = 2;
        for (let attempt = 1; attempt <= maxAttempts; attempt += 1) {
            try {
                const controller = new AbortController();
                const timeoutId = setTimeout(() => controller.abort(), 8000);
                const response = await fetch(url, {
                    method: 'GET',
                    headers: {
                        'Content-Type': 'application/json',
                    },
                    signal: controller.signal,
                });
                clearTimeout(timeoutId);
                console.log(`[API] Health check response: ${response.status} ${response.statusText}`);
                return response.ok;
            } catch (error) {
                console.error(`[API] Health check failed (attempt ${attempt}/${maxAttempts}):`, error.name, error.message);
                if (attempt === maxAttempts) {
                    return false;
                }
            }
        }
        return false;
    }

    /**
     * 计算散射过程
     * @param {Object} params - 计算参数
     * @param {number} params.p1x - 入射粒子1的x方向动量
     * @param {number} params.p1y - 入射粒子1的y方向动量
     * @param {number} params.p1z - 入射粒子1的z方向动量
     * @param {number} params.p2x - 入射粒子2的x方向动量
     * @param {number} params.p2y - 入射粒子2的y方向动量
     * @param {number} params.p2z - 入射粒子2的z方向动量
     * @param {number} params.m1 - 粒子1质量
     * @param {number} params.m2 - 粒子2质量
     * @param {number} params.m3 - 粒子3质量
     * @param {number} params.m4 - 粒子4质量
     * @param {number} [params.theta_star] - 质心系极角（可选）
     * @param {number} [params.phi_star] - 质心系方位角（可选）
     * @returns {Promise<Object>} 计算结果
     */
    static async computeScattering(params) {
        try {
            const response = await fetch(build_api_url('/compute'), {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(params),
            });

            if (!response.ok) {
                const errorData = await response.json();
                throw new Error(errorData.error || `HTTP error! status: ${response.status}`);
            }

            const result = await response.json();
            
            if (!result.success) {
                throw new Error(result.error || 'Computation failed');
            }

            return result.data;
        } catch (error) {
            console.error('API call failed:', error);
            throw error;
        }
    }

    static async createScanJob(payload) {
        const validation = validate_scan_request(payload);
        if (!validation.valid) {
            throw normalize_api_error({
                status: 400,
                payload: {
                    code: 'INVALID_REQUEST',
                    diagnostics: { errors: validation.errors },
                },
                message: validation.errors.join('\n'),
            });
        }

        return await request_json('/api/modules/pnjl-scan/jobs', {
            method: 'POST',
            body: JSON.stringify(payload),
            timeoutMs: 20000,
            retries: 0,
        });
    }

    static async getJobStatus(jobId) {
        if (!jobId) {
            throw normalize_api_error({
                status: 400,
                payload: { code: 'INVALID_REQUEST' },
                message: 'jobId 不能为空',
            });
        }

        return await request_json(`/api/modules/pnjl-scan/jobs/${encodeURIComponent(jobId)}`, {
            method: 'GET',
            timeoutMs: 10000,
            retries: 1,
        });
    }

    static async getJobResult(jobId) {
        if (!jobId) {
            throw normalize_api_error({
                status: 400,
                payload: { code: 'INVALID_REQUEST' },
                message: 'jobId 不能为空',
            });
        }

        return await request_json(`/api/modules/pnjl-scan/jobs/${encodeURIComponent(jobId)}/result`, {
            method: 'GET',
            timeoutMs: 10000,
            retries: 1,
        });
    }

    /**
     * 格式化错误消息
     * @param {Error} error - 错误对象
     * @returns {string} 格式化后的错误消息
     */
    static formatError(error) {
        if (error && error.code === 'QUEUE_FULL') {
            return '队列已满，请稍后重试';
        }
        if (error && error.code === 'INVALID_REQUEST') {
            return error.message || '请求参数无效，请检查输入';
        }
        if (error.name === 'TypeError' && error.message.includes('fetch')) {
            return '无法连接到服务器，请确保Julia服务器正在运行 (julia server.jl)';
        }
        if (error.name === 'AbortError') {
            return '请求超时，请检查服务器是否正常运行';
        }
        return error.message || '未知错误';
    }
}

/**
 * 验证输入参数
 * @param {Object} params - 输入参数
 * @returns {Object} { valid: boolean, errors: string[] }
 */
export function validateInput(params) {
    const errors = [];
    
    // 检查动量
    const momentumKeys = ['p1x', 'p1y', 'p1z', 'p2x', 'p2y', 'p2z'];
    for (const key of momentumKeys) {
        if (isNaN(params[key]) || !isFinite(params[key])) {
            errors.push(`${key} 必须是有效的数字`);
        }
    }
    
    // 检查质量（必须为正）
    const massKeys = ['m1', 'm2', 'm3', 'm4'];
    for (const key of massKeys) {
        if (isNaN(params[key]) || params[key] <= 0) {
            errors.push(`${key} 必须是正数`);
        }
    }
    
    // 检查角度范围
    if (params.theta_star !== undefined) {
        if (isNaN(params.theta_star) || params.theta_star < 0 || params.theta_star > Math.PI) {
            errors.push('theta_star 必须在 [0, π] 范围内');
        }
    }
    
    if (params.phi_star !== undefined) {
        if (isNaN(params.phi_star) || params.phi_star < 0 || params.phi_star > 2 * Math.PI) {
            errors.push('phi_star 必须在 [0, 2π] 范围内');
        }
    }
    
    return {
        valid: errors.length === 0,
        errors: errors
    };
}

/**
 * 格式化向量显示
 * @param {Array<number>} vector - 三维向量
 * @param {number} [precision=4] - 小数精度
 * @returns {string} 格式化字符串
 */
export function formatVector(vector, precision = 4) {
    if (!Array.isArray(vector) || vector.length !== 3) {
        return 'N/A';
    }
    return `[${vector.map(v => v.toFixed(precision)).join(', ')}]`;
}

/**
 * 格式化标量显示
 * @param {number} value - 标量值
 * @param {number} [precision=4] - 小数精度
 * @returns {string} 格式化字符串
 */
export function formatScalar(value, precision = 4) {
    if (isNaN(value) || !isFinite(value)) {
        return 'N/A';
    }
    return value.toFixed(precision);
}

/**
 * 检查数值是否接近零
 * @param {number} value - 数值
 * @param {number} [tolerance=1e-9] - 容差
 * @returns {boolean}
 */
export function isNearZero(value, tolerance = 1e-9) {
    return Math.abs(value) < tolerance;
}
