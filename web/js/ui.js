/**
 * UI控制模块
 * 负责用户交互和界面更新
 */

import { API, validateInput, formatVector, formatScalar, isNearZero } from './api.js';
import { get_api_base_url, set_api_base_url } from './runtime_config.js';
import { Visualization } from './visualization.js';

export class UI {
    constructor() {
        this.visualization = null;
        this.isComputing = false;
        this.pollTimer = null;
        this.scanRunning = false;
        this.activeJobId = null;
        this.lastScanPayload = null;
        
        // DOM元素
        this.form = null;
        this.computeBtn = null;
        this.resultsPanel = null;
        this.statusIndicator = null;
        this.statusText = null;

        this.scanForm = null;
        this.scanSubmitBtn = null;
        this.scanRetryBtn = null;
        this.scanStatus = null;
        this.scanDetails = null;
        this.scanResultPanel = null;
        this.scanTechError = null;
        this.scanTechDetails = null;
        this.apiBaseUrlInput = null;
        this.apiBaseUrlBtn = null;
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
        this.scanForm = document.getElementById('scan-form');
        this.scanSubmitBtn = document.getElementById('scan-submit-btn');
        this.scanRetryBtn = document.getElementById('scan-retry-btn');
        this.scanStatus = document.getElementById('scan-job-status');
        this.scanDetails = document.getElementById('scan-job-details');
        this.scanResultPanel = document.getElementById('scan-result-panel');
        this.scanTechError = document.getElementById('scan-tech-error');
        this.scanTechDetails = document.getElementById('scan-tech-details');
        this.apiBaseUrlInput = document.getElementById('api-base-url');
        this.apiBaseUrlBtn = document.getElementById('apply-api-base-url');

        // 初始化可视化
        this.visualization = new Visualization('canvas-container');

        if (this.apiBaseUrlInput) {
            this.apiBaseUrlInput.value = get_api_base_url();
        }
        
        // 绑定事件
        this.form.addEventListener('submit', (e) => this.handleSubmit(e));
        if (this.scanForm) {
            this.scanForm.addEventListener('submit', (e) => this.handleScanSubmit(e));
        }
        if (this.scanRetryBtn) {
            this.scanRetryBtn.addEventListener('click', () => this.handleScanRetry());
        }
        if (this.apiBaseUrlBtn) {
            this.apiBaseUrlBtn.addEventListener('click', () => this.applyApiBaseUrl());
        }
        const scanKind = document.getElementById('scan-kind');
        if (scanKind) {
            scanKind.addEventListener('change', () => this.syncScanKindRows());
            this.syncScanKindRows();
        }
        
        // 检查服务器状态
        this.checkServerStatus();
        
        // 定期检查服务器状态
        setInterval(() => this.checkServerStatus(), 10000);
    }

    syncScanKindRows() {
        const kind = this.getScanKind();
        const muRow = document.getElementById('scan-mu-row');
        const rhoRow = document.getElementById('scan-rho-row');
        if (!muRow || !rhoRow) {
            return;
        }
        if (kind === 'trho') {
            muRow.style.display = 'none';
            rhoRow.style.display = '';
        } else {
            muRow.style.display = '';
            rhoRow.style.display = 'none';
        }
    }

    getScanKind() {
        const kindEl = document.getElementById('scan-kind');
        return kindEl ? String(kindEl.value).toLowerCase() : 'tmu';
    }

    parseNumberList(raw, fieldName, required = true) {
        const text = String(raw || '').trim();
        if (!text) {
            if (required) {
                throw new Error(`${fieldName} 不能为空`);
            }
            return [];
        }
        const values = text.split(',').map((item) => Number(item.trim()));
        if (!values.every((v) => Number.isFinite(v))) {
            throw new Error(`${fieldName} 需要逗号分隔数字列表`);
        }
        return values;
    }

    collectScanPayload() {
        const kind = this.getScanKind();
        const tValues = this.parseNumberList(document.getElementById('scan-t-values')?.value, 'T_values');

        const payload = {
            kind,
            params: {
                T_values: tValues,
                max_retries: 0,
            },
        };

        if (kind === 'trho') {
            payload.params.rho_values = this.parseNumberList(document.getElementById('scan-rho-values')?.value, 'rho_values');
        } else {
            payload.params.mu_values = this.parseNumberList(document.getElementById('scan-mu-values')?.value, 'mu_values');
        }

        const selectedStrategy = document.querySelector('input[name="xi-strategy"]:checked')?.value;
        if (!selectedStrategy) {
            throw new Error('请先选择 xi 策略');
        }

        if (selectedStrategy === 'xi') {
            payload.params.xi = Number(document.getElementById('scan-xi')?.value);
        } else if (selectedStrategy === 'xi_values') {
            payload.params.xi_values = this.parseNumberList(document.getElementById('scan-xi-values')?.value, 'xi_values');
        } else if (selectedStrategy === 'xi_grid') {
            const start = Number(document.getElementById('scan-xi-grid-start')?.value);
            const stop = Number(document.getElementById('scan-xi-grid-stop')?.value);
            const step = Number(document.getElementById('scan-xi-grid-step')?.value);
            payload.params.xi_grid = { start, stop, step };
        }

        return payload;
    }

    applyApiBaseUrl() {
        if (!this.apiBaseUrlInput) {
            return;
        }
        const updated = set_api_base_url(this.apiBaseUrlInput.value);
        this.apiBaseUrlInput.value = updated;
        this.showSuccess(`API_BASE_URL 已更新为 ${updated}`);
        this.checkServerStatus();
    }

    setScanStatus(text, level = 'info') {
        if (!this.scanStatus) {
            return;
        }
        this.scanStatus.textContent = text;
        this.scanStatus.className = `scan-status ${level}`;
    }

    writeScanDetails(data) {
        if (!this.scanDetails) {
            return;
        }
        this.scanDetails.textContent = typeof data === 'string' ? data : JSON.stringify(data, null, 2);
    }

    clearTechError() {
        if (this.scanTechError) {
            this.scanTechError.style.display = 'none';
        }
        if (this.scanTechDetails) {
            this.scanTechDetails.textContent = '-';
        }
    }

    showTechError(error) {
        if (!this.scanTechError || !this.scanTechDetails) {
            return;
        }
        this.scanTechError.style.display = '';
        this.scanTechDetails.textContent = JSON.stringify({
            code: error?.code || null,
            status: error?.status || null,
            message: error?.message || 'unknown',
            diagnostics: error?.diagnostics || null,
            payload: error?.payload || null,
        }, null, 2);
    }

    updateScanResult(resultPayload) {
        const stats = resultPayload?.result?.stats || {};
        const outputPath = resultPayload?.result?.output_path || '-';
        document.getElementById('scan-total').textContent = stats.total ?? '-';
        document.getElementById('scan-success').textContent = stats.success ?? '-';
        document.getElementById('scan-failure').textContent = stats.failure ?? '-';
        document.getElementById('scan-skipped').textContent = stats.skipped ?? '-';
        document.getElementById('scan-output-path').textContent = outputPath;
        if (this.scanResultPanel) {
            this.scanResultPanel.style.display = '';
        }
    }

    stopPolling() {
        if (this.pollTimer) {
            clearInterval(this.pollTimer);
            this.pollTimer = null;
        }
        this.scanRunning = false;
        if (this.scanSubmitBtn) {
            this.scanSubmitBtn.disabled = false;
            this.scanSubmitBtn.textContent = '提交扫描任务';
        }
    }

    startPolling(jobId) {
        this.stopPolling();
        this.scanRunning = true;
        if (this.scanSubmitBtn) {
            this.scanSubmitBtn.disabled = true;
            this.scanSubmitBtn.textContent = '轮询中...';
        }

        this.pollTimer = setInterval(async () => {
            try {
                const statusPayload = await API.getJobStatus(jobId);
                this.writeScanDetails(statusPayload);
                const jobStatus = statusPayload?.job_status;
                const percent = statusPayload?.progress?.percent;
                if (jobStatus === 'queued') {
                    this.setScanStatus(`任务排队中 (${percent ?? 0}%)`, 'queued');
                } else if (jobStatus === 'running') {
                    this.setScanStatus(`任务运行中 (${percent ?? 0}%)`, 'running');
                } else if (jobStatus === 'succeeded') {
                    this.setScanStatus('任务成功，正在加载结果...', 'success');
                    this.stopPolling();
                    const resultPayload = await API.getJobResult(jobId);
                    this.updateScanResult(resultPayload);
                    this.writeScanDetails(resultPayload);
                    this.clearTechError();
                    if (this.scanRetryBtn) {
                        this.scanRetryBtn.disabled = true;
                    }
                } else if (jobStatus === 'failed') {
                    this.stopPolling();
                    this.setScanStatus('任务失败，可查看详情后重试', 'error');
                    if (this.scanRetryBtn) {
                        this.scanRetryBtn.disabled = false;
                    }
                }
            } catch (error) {
                this.stopPolling();
                this.setScanStatus(API.formatError(error), 'error');
                this.showTechError(error);
                if (this.scanRetryBtn) {
                    this.scanRetryBtn.disabled = false;
                }
            }
        }, 3000);
    }

    async handleScanSubmit(event) {
        event.preventDefault();
        if (this.scanRunning) {
            return;
        }

        this.clearError();
        this.clearTechError();
        if (this.scanResultPanel) {
            this.scanResultPanel.style.display = 'none';
        }
        try {
            const payload = this.collectScanPayload();
            this.lastScanPayload = payload;
            const created = await API.createScanJob(payload);
            this.activeJobId = created.job_id;
            this.setScanStatus(`任务已创建: ${created.job_id}`, 'queued');
            this.writeScanDetails(created);
            this.startPolling(created.job_id);
        } catch (error) {
            this.setScanStatus(API.formatError(error), 'error');
            this.showError(API.formatError(error));
            this.showTechError(error);
            if (this.scanRetryBtn) {
                this.scanRetryBtn.disabled = !this.lastScanPayload;
            }
        }
    }

    async handleScanRetry() {
        if (!this.lastScanPayload || this.scanRunning) {
            return;
        }
        this.clearTechError();
        try {
            const created = await API.createScanJob(this.lastScanPayload);
            this.activeJobId = created.job_id;
            this.setScanStatus(`重试已创建任务: ${created.job_id}`, 'queued');
            this.writeScanDetails(created);
            this.startPolling(created.job_id);
        } catch (error) {
            this.setScanStatus(API.formatError(error), 'error');
            this.showTechError(error);
        }
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
            this.statusText.textContent = `服务器在线 (${get_api_base_url()})`;
            this.computeBtn.disabled = false;
        } else {
            this.statusIndicator.className = 'status-indicator offline';
            this.statusText.textContent = `服务器离线 (${get_api_base_url()}) - 请运行 julia scripts/server/server_full.jl`;
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
