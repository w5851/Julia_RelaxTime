/**
 * UI控制模块
 * 负责用户交互和界面更新
 */

import { API, validateInput, formatVector, formatScalar, isNearZero } from './api.js';
import { Visualization } from './visualization.js';

export class UI {
    constructor() {
        this.visualization = null;
        this.isComputing = false;
        
        // DOM元素
        this.form = null;
        this.computeBtn = null;
        this.resultsPanel = null;
        this.statusIndicator = null;
        this.statusText = null;
    }

    /**
     * 初始化UI
     */
    initialize() {
        // 获取DOM元素
        this.form = document.getElementById('input-form');
        this.computeBtn = document.getElementById('compute-btn');
        this.resultsPanel = document.getElementById('results-panel');
        this.statusIndicator = document.getElementById('server-status');
        this.statusText = document.getElementById('status-text');
        
        // 初始化可视化
        this.visualization = new Visualization('canvas-container');
        
        // 绑定事件
        this.form.addEventListener('submit', (e) => this.handleSubmit(e));
        
        // 检查服务器状态
        this.checkServerStatus();
        
        // 定期检查服务器状态
        setInterval(() => this.checkServerStatus(), 10000);
    }

    /**
     * 检查服务器状态
     */
    async checkServerStatus() {
        this.statusIndicator.className = 'status-indicator checking';
        this.statusText.textContent = '检查服务器...';
        
        console.log('[UI] Checking server health...');
        const isOnline = await API.checkHealth();
        console.log('[UI] Server status:', isOnline ? 'ONLINE' : 'OFFLINE');
        
        if (isOnline) {
            this.statusIndicator.className = 'status-indicator online';
            this.statusText.textContent = '服务器在线';
            this.computeBtn.disabled = false;
        } else {
            this.statusIndicator.className = 'status-indicator offline';
            this.statusText.textContent = '服务器离线 - 请运行 julia server.jl';
            this.computeBtn.disabled = true;
        }
    }

    /**
     * 处理表单提交
     * @param {Event} event - 提交事件
     */
    async handleSubmit(event) {
        event.preventDefault();
        
        if (this.isComputing) {
            return;
        }
        
        // 收集输入参数
        const params = this.collectFormData();
        
        // 验证输入
        const validation = validateInput(params);
        if (!validation.valid) {
            this.showError(validation.errors.join('\n'));
            return;
        }
        
        // 开始计算
        this.isComputing = true;
        this.computeBtn.classList.add('loading');
        this.computeBtn.textContent = '计算中...';
        this.clearError();
        
        try {
            // 调用API
            const result = await API.computeScattering(params);
            
            // 更新可视化
            this.visualization.update(result);
            
            // 更新结果显示
            this.updateResults(result);
            
            // 显示结果面板
            this.resultsPanel.style.display = 'block';
            
        } catch (error) {
            console.error('Computation failed:', error);
            this.showError(API.formatError(error));
        } finally {
            this.isComputing = false;
            this.computeBtn.classList.remove('loading');
            this.computeBtn.textContent = '🚀 计算散射';
        }
    }

    /**
     * 收集表单数据
     * @returns {Object} 表单参数
     */
    collectFormData() {
        const getNumber = (id) => {
            const element = document.getElementById(id);
            const value = parseFloat(element.value);
            console.log(`[UI] Collected ${id} = ${value} (raw: "${element.value}")`);
            return value;
        };
        
        const params = {
            p1x: getNumber('p1x'),
            p1y: getNumber('p1y'),
            p1z: getNumber('p1z'),
            p2x: getNumber('p2x'),
            p2y: getNumber('p2y'),
            p2z: getNumber('p2z'),
            m1: getNumber('m1'),
            m2: getNumber('m2'),
            m3: getNumber('m3'),
            m4: getNumber('m4'),
            theta_star: getNumber('theta_star'),
            phi_star: getNumber('phi_star')
        };
        
        console.log('[UI] Collected form data:', params);
        return params;
    }

    /**
     * 更新结果显示
     * @param {Object} data - 计算结果
     */
    updateResults(data) {
        const { physics, momenta, ellipsoid, validation } = data;
        
        // 更新物理量
        document.getElementById('sqrt-s').textContent = formatScalar(physics.sqrt_s);
        document.getElementById('p-star').textContent = formatScalar(physics.p_star);
        document.getElementById('beta').textContent = formatScalar(physics.beta);
        document.getElementById('gamma').textContent = formatScalar(physics.gamma);
        
        // 更新出射动量
        document.getElementById('p3-lab').textContent = formatVector(momenta.p3);
        document.getElementById('p4-lab').textContent = formatVector(momenta.p4);
        
        // 更新椭球参数
        document.getElementById('ellipsoid-center').textContent = formatVector(ellipsoid.center);
        document.getElementById('ellipsoid-lengths').textContent = formatVector(ellipsoid.half_lengths);
        
        // 更新验证结果
        this.updateValidation(validation);
    }

    /**
     * 更新物理验证显示
     * @param {Object} validation - 验证结果
     */
    updateValidation(validation) {
        const energyCheck = document.getElementById('energy-check');
        const momentumCheck = document.getElementById('momentum-check');
        
        const energyValid = isNearZero(validation.energy_conservation);
        const momentumValid = isNearZero(validation.momentum_conservation);
        
        // 能量守恒
        energyCheck.textContent = energyValid ? '✓ 通过' : `✗ 失败 (ΔE=${formatScalar(validation.energy_conservation, 6)})`;
        energyCheck.className = energyValid ? 'value valid' : 'value invalid';
        
        // 动量守恒
        momentumCheck.textContent = momentumValid ? '✓ 通过' : `✗ 失败 (|Δp|=${formatScalar(validation.momentum_conservation, 6)})`;
        momentumCheck.className = momentumValid ? 'value valid' : 'value invalid';
    }

    /**
     * 显示错误消息
     * @param {string} message - 错误消息
     */
    showError(message) {
        // 移除旧错误消息
        this.clearError();
        
        // 创建错误元素
        const errorDiv = document.createElement('div');
        errorDiv.className = 'error-message';
        errorDiv.id = 'error-message';
        errorDiv.textContent = message;
        
        // 插入到表单后
        this.form.insertAdjacentElement('afterend', errorDiv);
        
        // 5秒后自动消失
        setTimeout(() => this.clearError(), 5000);
    }

    /**
     * 清除错误消息
     */
    clearError() {
        const errorDiv = document.getElementById('error-message');
        if (errorDiv) {
            errorDiv.remove();
        }
    }

    /**
     * 显示成功消息
     * @param {string} message - 成功消息
     */
    showSuccess(message) {
        const successDiv = document.createElement('div');
        successDiv.className = 'success-message';
        successDiv.textContent = message;
        
        this.form.insertAdjacentElement('afterend', successDiv);
        
        setTimeout(() => successDiv.remove(), 3000);
    }
}
