/**
 * API通信模块
 * 负责与Julia HTTP服务器通信
 */

const API_BASE_URL = 'http://localhost:8080';

export class API {
    /**
     * 检查服务器健康状态
     * @returns {Promise<boolean>} 服务器是否在线
     */
    static async checkHealth() {
        try {
            console.log(`[API] Checking health at ${API_BASE_URL}/health`);
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 3000); // 3秒超时
            
            const response = await fetch(`${API_BASE_URL}/health`, {
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
            console.error('[API] Health check failed:', error.name, error.message);
            return false;
        }
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
            const controller = new AbortController();
            const timeoutId = setTimeout(() => controller.abort(), 10000); // 10秒超时
            
            const response = await fetch(`${API_BASE_URL}/compute`, {
                method: 'POST',
                headers: {
                    'Content-Type': 'application/json',
                },
                body: JSON.stringify(params),
                signal: controller.signal,
            });
            
            clearTimeout(timeoutId);

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

    /**
     * 格式化错误消息
     * @param {Error} error - 错误对象
     * @returns {string} 格式化后的错误消息
     */
    static formatError(error) {
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
