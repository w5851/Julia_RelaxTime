/**
 * UI控制模块
 * 负责用户交互和界面更新
 */

import { API, validateInput, formatVector, formatScalar, isNearZero } from './api.js';
import { get_api_base_url, set_api_base_url } from './runtime_config.js';
import { Visualization } from './visualization.js';
import { get_default_templates, load_task_history, load_task_templates, upsert_task_history_entry } from './task_store.js';

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
        this.scanCancelBtn = null;
        this.scanTemplateSelect = null;
        this.scanTemplateApplyBtn = null;
        this.scanJobIdInput = null;
        this.scanJobLoadBtn = null;
        this.scanHistoryList = null;
        this.scanResultPreview = null;
        this.downloadResultBtn = null;
        this.copyParamsBtn = null;
        this.copyCurlBtn = null;
        this.scanReproPanel = null;
        this.scanReproPreview = null;
        this.lastResultPayload = null;
        this.progressChart = null;
        this.historyHeatmap = null;
        this.navScatteringBtn = null;
        this.navTaskCenterBtn = null;
        this.serviceTabs = [];
        this.servicePanels = [];
        this.pnjlGapForm = null;
        this.pnjlGapSubmitBtn = null;
        this.pnjlGapCopyParamsBtn = null;
        this.pnjlGapCopyCurlBtn = null;
        this.pnjlGapStatus = null;
        this.pnjlGapResultPanel = null;
        this.pnjlGapResultPreview = null;
        this.lastPnjlGapPayload = null;
        this.lastPnjlGapResultPayload = null;
        this.transportPointForm = null;
        this.transportPointSubmitBtn = null;
        this.transportPointCopyParamsBtn = null;
        this.transportPointCopyCurlBtn = null;
        this.transportPointStatus = null;
        this.transportPointResultPanel = null;
        this.transportPointResultPreview = null;
        this.lastTransportPointPayload = null;
        this.lastTransportPointResultPayload = null;
        this.scriptTaskCatalog = [];
        this.activeScriptTaskId = null;
        this.activeScriptTaskJobId = null;
        this.scriptTaskPollTimer = null;
        this.scriptTaskRunning = false;
        this.lastScriptTaskPayload = null;
        this.lastScriptTaskResultPayload = null;
        this.scriptTaskCatalogEl = null;
        this.scriptTaskRefreshBtn = null;
        this.scriptTaskForm = null;
        this.scriptTaskSelect = null;
        this.scriptTaskPresetSelect = null;
        this.scriptTaskConfirmHeavy = null;
        this.scriptTaskCustomFieldset = null;
        this.scriptTaskCustomArgs = null;
        this.scriptTaskSelectedSummary = null;
        this.scriptTaskSubmitBtn = null;
        this.scriptTaskCancelBtn = null;
        this.scriptTaskCopyCurlBtn = null;
        this.scriptTaskStatus = null;
        this.scriptTaskDetails = null;
        this.scriptTaskResultPanel = null;
        this.scriptTaskResultPreview = null;
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
        this.scanCancelBtn = document.getElementById('scan-cancel-btn');
        this.scanTemplateSelect = document.getElementById('scan-template-select');
        this.scanTemplateApplyBtn = document.getElementById('scan-template-apply-btn');
        this.scanJobIdInput = document.getElementById('scan-job-id-input');
        this.scanJobLoadBtn = document.getElementById('scan-job-load-btn');
        this.scanHistoryList = document.getElementById('scan-history-list');
        this.scanResultPreview = document.getElementById('scan-result-preview');
        this.downloadResultBtn = document.getElementById('scan-download-btn');
        this.copyParamsBtn = document.getElementById('scan-copy-params-btn');
        this.copyCurlBtn = document.getElementById('scan-copy-curl-btn');
        this.scanReproPanel = document.getElementById('scan-repro-panel');
        this.scanReproPreview = document.getElementById('scan-repro-preview');
        this.progressChart = document.getElementById('scan-progress-chart');
        this.historyHeatmap = document.getElementById('scan-history-heatmap');
        this.navScatteringBtn = document.getElementById('nav-scattering');
        this.navTaskCenterBtn = document.getElementById('nav-task-center');
        this.serviceTabs = Array.from(document.querySelectorAll('[data-service-tab]'));
        this.servicePanels = Array.from(document.querySelectorAll('[data-service-panel]'));
        this.pnjlGapForm = document.getElementById('pnjl-gap-form');
        this.pnjlGapSubmitBtn = document.getElementById('pnjl-gap-submit-btn');
        this.pnjlGapCopyParamsBtn = document.getElementById('pnjl-gap-copy-params-btn');
        this.pnjlGapCopyCurlBtn = document.getElementById('pnjl-gap-copy-curl-btn');
        this.pnjlGapStatus = document.getElementById('pnjl-gap-status');
        this.pnjlGapResultPanel = document.getElementById('pnjl-gap-result-panel');
        this.pnjlGapResultPreview = document.getElementById('pnjl-gap-result-preview');
        this.transportPointForm = document.getElementById('transport-point-form');
        this.transportPointSubmitBtn = document.getElementById('transport-point-submit-btn');
        this.transportPointCopyParamsBtn = document.getElementById('transport-point-copy-params-btn');
        this.transportPointCopyCurlBtn = document.getElementById('transport-point-copy-curl-btn');
        this.transportPointStatus = document.getElementById('transport-point-status');
        this.transportPointResultPanel = document.getElementById('transport-point-result-panel');
        this.transportPointResultPreview = document.getElementById('transport-point-result-preview');
        this.scriptTaskCatalogEl = document.getElementById('script-task-catalog');
        this.scriptTaskRefreshBtn = document.getElementById('script-task-refresh-btn');
        this.scriptTaskForm = document.getElementById('script-task-form');
        this.scriptTaskSelect = document.getElementById('script-task-select');
        this.scriptTaskPresetSelect = document.getElementById('script-task-preset');
        this.scriptTaskConfirmHeavy = document.getElementById('script-task-confirm-heavy');
        this.scriptTaskCustomFieldset = document.getElementById('script-task-custom-fieldset');
        this.scriptTaskCustomArgs = document.getElementById('script-task-custom-args');
        this.scriptTaskSelectedSummary = document.getElementById('script-task-selected-summary');
        this.scriptTaskSubmitBtn = document.getElementById('script-task-submit-btn');
        this.scriptTaskCancelBtn = document.getElementById('script-task-cancel-btn');
        this.scriptTaskCopyCurlBtn = document.getElementById('script-task-copy-curl-btn');
        this.scriptTaskStatus = document.getElementById('script-task-status');
        this.scriptTaskDetails = document.getElementById('script-task-details');
        this.scriptTaskResultPanel = document.getElementById('script-task-result-panel');
        this.scriptTaskResultPreview = document.getElementById('script-task-result-preview');

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
        if (this.scanCancelBtn) {
            this.scanCancelBtn.addEventListener('click', () => this.handleScanCancel());
        }
        if (this.scanTemplateApplyBtn) {
            this.scanTemplateApplyBtn.addEventListener('click', () => this.applyTemplateSelection());
        }
        if (this.scanJobLoadBtn) {
            this.scanJobLoadBtn.addEventListener('click', () => this.handleLoadTaskByJobId());
        }
        if (this.downloadResultBtn) {
            this.downloadResultBtn.addEventListener('click', () => this.handleDownloadResult());
        }
        if (this.copyParamsBtn) {
            this.copyParamsBtn.addEventListener('click', () => this.handleCopyScanParams());
        }
        if (this.copyCurlBtn) {
            this.copyCurlBtn.addEventListener('click', () => this.handleCopyScanCurl());
        }
        if (this.navScatteringBtn) {
            this.navScatteringBtn.addEventListener('click', () => this.navigateToPage('scattering'));
        }
        if (this.navTaskCenterBtn) {
            this.navTaskCenterBtn.addEventListener('click', () => this.navigateToPage('task-center'));
        }
        this.serviceTabs.forEach((tab) => {
            tab.addEventListener('click', () => this.switchServicePanel(tab.dataset.serviceTab));
        });
        if (this.pnjlGapForm) {
            this.pnjlGapForm.addEventListener('submit', (e) => this.handlePnjlGapSubmit(e));
        }
        if (this.pnjlGapCopyParamsBtn) {
            this.pnjlGapCopyParamsBtn.addEventListener('click', () => this.handleCopyPnjlGapParams());
        }
        if (this.pnjlGapCopyCurlBtn) {
            this.pnjlGapCopyCurlBtn.addEventListener('click', () => this.handleCopyPnjlGapCurl());
        }
        if (this.transportPointForm) {
            this.transportPointForm.addEventListener('submit', (e) => this.handleTransportPointSubmit(e));
        }
        if (this.transportPointCopyParamsBtn) {
            this.transportPointCopyParamsBtn.addEventListener('click', () => this.handleCopyTransportPointParams());
        }
        if (this.transportPointCopyCurlBtn) {
            this.transportPointCopyCurlBtn.addEventListener('click', () => this.handleCopyTransportPointCurl());
        }
        if (this.scriptTaskRefreshBtn) {
            this.scriptTaskRefreshBtn.addEventListener('click', () => this.initializeScriptTasks());
        }
        if (this.scriptTaskForm) {
            this.scriptTaskForm.addEventListener('submit', (e) => this.handleScriptTaskSubmit(e));
        }
        if (this.scriptTaskSelect) {
            this.scriptTaskSelect.addEventListener('change', () => this.selectScriptTask(this.scriptTaskSelect.value));
        }
        if (this.scriptTaskPresetSelect) {
            this.scriptTaskPresetSelect.addEventListener('change', () => this.syncScriptTaskPresetState());
        }
        if (this.scriptTaskCancelBtn) {
            this.scriptTaskCancelBtn.addEventListener('click', () => this.handleScriptTaskCancel());
        }
        if (this.scriptTaskCopyCurlBtn) {
            this.scriptTaskCopyCurlBtn.addEventListener('click', () => this.handleCopyScriptTaskCurl());
        }
        if (typeof window !== 'undefined') {
            window.addEventListener('hashchange', () => this.applyNavigationRoute());
        }
        const scanKind = document.getElementById('scan-kind');
        const scanMode = document.getElementById('scan-mode');
        if (scanKind) {
            scanKind.addEventListener('change', () => this.syncScanKindRows());
        }
        if (scanMode) {
            scanMode.addEventListener('change', () => this.syncScanFormRows());
        }
        this.syncScanFormRows();

        this.initializeTemplates();
        this.initializeScriptTasks();
        this.renderTaskHistory();
        this.applyNavigationRoute();
        
        // 检查服务器状态
        this.checkServerStatus();
        
        // 定期检查服务器状态
        setInterval(() => this.checkServerStatus(), 10000);
    }

    syncScanKindRows() {
        const kind = this.getScanKind();
        const muRow = document.getElementById('scan-mu-row');
        const rhoRow = document.getElementById('scan-rho-row');
        const pointMuRow = document.getElementById('scan-point-mu-row');
        const pointRhoRow = document.getElementById('scan-point-rho-row');
        if (!muRow || !rhoRow) {
            return;
        }
        if (kind === 'trho') {
            muRow.style.display = 'none';
            rhoRow.style.display = '';
            if (pointMuRow) {
                pointMuRow.style.display = 'none';
            }
            if (pointRhoRow) {
                pointRhoRow.style.display = '';
            }
        } else {
            muRow.style.display = '';
            rhoRow.style.display = 'none';
            if (pointMuRow) {
                pointMuRow.style.display = '';
            }
            if (pointRhoRow) {
                pointRhoRow.style.display = 'none';
            }
        }
    }

    getScanMode() {
        const modeEl = document.getElementById('scan-mode');
        return modeEl ? String(modeEl.value).toLowerCase() : 'scan';
    }

    syncScanFormRows() {
        const mode = this.getScanMode();
        const gridFieldset = document.getElementById('scan-grid-fieldset');
        const pointFieldset = document.getElementById('scan-point-fieldset');
        if (gridFieldset) {
            gridFieldset.style.display = mode === 'point' ? 'none' : '';
        }
        if (pointFieldset) {
            pointFieldset.style.display = mode === 'point' ? '' : 'none';
        }
        this.syncScanKindRows();
    }

    getScanKind() {
        const kindEl = document.getElementById('scan-kind');
        return kindEl ? String(kindEl.value).toLowerCase() : 'tmu';
    }

    switchServicePanel(panelId) {
        const target = panelId || 'scan';
        this.serviceTabs.forEach((tab) => {
            tab.classList.toggle('active', tab.dataset.serviceTab === target);
        });
        this.servicePanels.forEach((panel) => {
            panel.classList.toggle('active', panel.dataset.servicePanel === target);
        });
    }

    async initializeScriptTasks() {
        if (!this.scriptTaskCatalogEl) {
            return;
        }
        this.scriptTaskCatalogEl.textContent = '正在加载任务目录...';
        try {
            const payload = await API.getScriptTaskCatalog();
            this.scriptTaskCatalog = Array.isArray(payload?.tasks) ? payload.tasks : [];
            this.renderScriptTaskCatalog();
            this.populateScriptTaskSelect();
            if (this.scriptTaskCatalog.length > 0) {
                this.selectScriptTask(this.activeScriptTaskId || this.scriptTaskCatalog[0].id);
            }
        } catch (error) {
            this.scriptTaskCatalogEl.textContent = API.formatError(error);
            this.setScriptTaskStatus(API.formatError(error), 'error');
        }
    }

    getScriptTaskById(taskId = this.activeScriptTaskId) {
        return this.scriptTaskCatalog.find((task) => String(task.id) === String(taskId)) || null;
    }

    renderScriptTaskCatalog() {
        if (!this.scriptTaskCatalogEl) {
            return;
        }
        this.scriptTaskCatalogEl.innerHTML = '';
        if (!this.scriptTaskCatalog.length) {
            this.scriptTaskCatalogEl.textContent = '暂无脚本任务';
            return;
        }

        this.scriptTaskCatalog.forEach((task) => {
            const card = document.createElement('button');
            card.type = 'button';
            card.className = 'script-task-card';
            card.dataset.taskId = task.id;

            const title = document.createElement('h3');
            title.textContent = task.name || task.id;
            card.appendChild(title);

            const script = document.createElement('code');
            script.textContent = task.script || '-';
            card.appendChild(script);

            const purpose = document.createElement('p');
            purpose.textContent = task.purpose || '未提供用途说明';
            card.appendChild(purpose);

            const addField = (label, value) => {
                const field = document.createElement('div');
                field.className = 'script-task-field';
                field.textContent = `${label}: ${value || '-'}`;
                card.appendChild(field);
            };

            const useCases = Array.isArray(task.use_cases) ? task.use_cases.slice(0, 3) : [];
            if (useCases.length > 0) {
                const list = document.createElement('ul');
                useCases.forEach((item) => {
                    const li = document.createElement('li');
                    li.textContent = item;
                    list.appendChild(li);
                });
                card.appendChild(list);
            }

            addField('关键参数', Array.isArray(task.key_params) ? task.key_params.join(', ') : '-');
            addField('输出产物', Array.isArray(task.outputs) ? task.outputs.join(', ') : '-');
            addField('预计耗时', this.formatScriptTaskEstimatedTime(task));
            addField('本机建议', task.local_recommendation || '-');

            const meta = document.createElement('div');
            meta.className = 'script-task-meta';
            const presetNames = Object.keys(task.presets || {}).join(' / ');
            meta.textContent = `默认 ${task.default_preset || 'smoke'} | presets: ${presetNames || '-'}`;
            card.appendChild(meta);

            card.addEventListener('click', () => this.selectScriptTask(task.id));
            this.scriptTaskCatalogEl.appendChild(card);
        });
    }

    populateScriptTaskSelect() {
        if (!this.scriptTaskSelect) {
            return;
        }
        this.scriptTaskSelect.innerHTML = '';
        this.scriptTaskCatalog.forEach((task) => {
            const option = document.createElement('option');
            option.value = task.id;
            option.textContent = `${task.id} - ${task.name || task.id}`;
            this.scriptTaskSelect.appendChild(option);
        });
    }

    selectScriptTask(taskId) {
        const task = this.getScriptTaskById(taskId);
        if (!task) {
            return;
        }
        this.activeScriptTaskId = task.id;
        if (this.scriptTaskSelect) {
            this.scriptTaskSelect.value = task.id;
        }
        if (this.scriptTaskCatalogEl) {
            Array.from(this.scriptTaskCatalogEl.querySelectorAll('.script-task-card')).forEach((card) => {
                card.classList.toggle('active', card.dataset.taskId === String(task.id));
            });
        }
        if (this.scriptTaskPresetSelect) {
            this.scriptTaskPresetSelect.innerHTML = '';
            Object.keys(task.presets || {}).forEach((preset) => {
                const option = document.createElement('option');
                option.value = preset;
                const heavy = task.presets[preset]?.heavy ? '重任务' : '安全预览';
                option.textContent = `${preset} (${heavy})`;
                this.scriptTaskPresetSelect.appendChild(option);
            });
            this.scriptTaskPresetSelect.value = task.default_preset || Object.keys(task.presets || {})[0] || 'smoke';
        }
        this.renderScriptTaskSummary(task);
        this.syncScriptTaskPresetState();
    }

    renderScriptTaskSummary(task = this.getScriptTaskById()) {
        if (!this.scriptTaskSelectedSummary || !task) {
            return;
        }
        const keyParams = Array.isArray(task.key_params) ? task.key_params.join(', ') : '-';
        const outputs = Array.isArray(task.outputs) ? task.outputs.join(', ') : '-';
        const estimatedTime = this.formatScriptTaskEstimatedTime(task);
        const recommendation = task.local_recommendation || '-';
        this.scriptTaskSelectedSummary.textContent = [
            `用途: ${task.purpose || '-'}`,
            `关键参数: ${keyParams}`,
            `输出: ${outputs}`,
            `预计耗时: ${estimatedTime}`,
            `本机建议: ${recommendation}`,
        ].join('\n');
    }

    formatScriptTaskEstimatedTime(task) {
        const estimated = task?.estimated_time || {};
        if (!estimated || typeof estimated !== 'object') {
            return '-';
        }
        const order = ['smoke', 'canonical', 'custom'];
        const keys = [
            ...order.filter((key) => estimated[key]),
            ...Object.keys(estimated).filter((key) => !order.includes(key)),
        ];
        if (keys.length === 0) {
            return '-';
        }
        return keys.map((key) => `${key}: ${estimated[key]}`).join('；');
    }

    syncScriptTaskPresetState() {
        const task = this.getScriptTaskById();
        const preset = String(this.scriptTaskPresetSelect?.value || task?.default_preset || 'smoke');
        const presetData = task?.presets?.[preset] || {};
        const isCustom = preset === 'custom';
        const isHeavy = !!presetData.heavy;
        if (this.scriptTaskCustomFieldset) {
            this.scriptTaskCustomFieldset.style.display = isCustom ? '' : 'none';
        }
        if (this.scriptTaskConfirmHeavy) {
            this.scriptTaskConfirmHeavy.disabled = !isHeavy;
            if (!isHeavy) {
                this.scriptTaskConfirmHeavy.checked = false;
            }
        }
    }

    readScriptTaskCustomArgs() {
        const text = String(this.scriptTaskCustomArgs?.value || '').trim();
        if (!text) {
            return [];
        }
        return text.split(/\r?\n/).map((line) => line.trim()).filter(Boolean);
    }

    collectScriptTaskPayload() {
        const task = this.getScriptTaskById();
        if (!task) {
            throw new Error('请先选择脚本任务');
        }
        const preset = String(this.scriptTaskPresetSelect?.value || task.default_preset || 'smoke');
        const payload = {
            task_id: task.id,
            preset,
            confirm_heavy: !!this.scriptTaskConfirmHeavy?.checked,
        };
        if (preset === 'custom') {
            payload.custom_args = this.readScriptTaskCustomArgs();
        }
        return payload;
    }

    buildScriptTaskCurlCommand() {
        if (!this.lastScriptTaskPayload) {
            return '';
        }
        const body = JSON.stringify(this.lastScriptTaskPayload);
        const escapedBody = body.replace(/'/g, "'\"'\"'");
        return `curl -X POST "${API.buildScriptTaskCreateUrl()}" -H "Content-Type: application/json" --data-raw '${escapedBody}'`;
    }

    setScriptTaskStatus(text, level = 'info') {
        if (!this.scriptTaskStatus) {
            return;
        }
        this.scriptTaskStatus.textContent = text;
        this.scriptTaskStatus.className = `scan-status ${level}`;
    }

    writeScriptTaskDetails(data) {
        if (!this.scriptTaskDetails) {
            return;
        }
        this.scriptTaskDetails.textContent = typeof data === 'string' ? data : JSON.stringify(data, null, 2);
    }

    stopScriptTaskPolling() {
        if (this.scriptTaskPollTimer) {
            clearInterval(this.scriptTaskPollTimer);
            this.scriptTaskPollTimer = null;
        }
        this.scriptTaskRunning = false;
        if (this.scriptTaskSubmitBtn) {
            this.scriptTaskSubmitBtn.disabled = false;
            this.scriptTaskSubmitBtn.textContent = '提交脚本任务';
        }
        if (this.scriptTaskCancelBtn) {
            this.scriptTaskCancelBtn.disabled = true;
        }
    }

    startScriptTaskPolling(jobId) {
        this.stopScriptTaskPolling();
        this.scriptTaskRunning = true;
        if (this.scriptTaskSubmitBtn) {
            this.scriptTaskSubmitBtn.disabled = true;
            this.scriptTaskSubmitBtn.textContent = '轮询中...';
        }
        if (this.scriptTaskCancelBtn) {
            this.scriptTaskCancelBtn.disabled = false;
        }

        this.scriptTaskPollTimer = setInterval(async () => {
            try {
                const statusPayload = await API.getScriptTaskJobStatus(jobId);
                this.writeScriptTaskDetails(statusPayload);
                const jobStatus = statusPayload?.job_status;
                if (jobStatus === 'queued') {
                    this.setScriptTaskStatus('任务排队中', 'queued');
                } else if (jobStatus === 'running') {
                    this.setScriptTaskStatus('任务运行中', 'running');
                } else if (['succeeded', 'failed', 'cancelled'].includes(jobStatus)) {
                    this.stopScriptTaskPolling();
                    const resultPayload = await API.getScriptTaskJobResult(jobId);
                    this.updateScriptTaskResult(resultPayload);
                    this.writeScriptTaskDetails(resultPayload);
                    const level = jobStatus === 'succeeded' ? 'success' : (jobStatus === 'failed' ? 'error' : 'info');
                    this.setScriptTaskStatus(`任务结束: ${jobStatus}`, level);
                }

                upsert_task_history_entry({
                    job_id: statusPayload.job_id,
                    kind: statusPayload.kind,
                    module: 'script-tasks',
                    task_id: statusPayload.task_id,
                    job_status: statusPayload.job_status,
                });
                this.renderTaskHistory();
            } catch (error) {
                this.stopScriptTaskPolling();
                this.setScriptTaskStatus(API.formatError(error), 'error');
                this.showTechError(error);
            }
        }, 3000);
    }

    updateScriptTaskResult(resultPayload) {
        const result = resultPayload?.result || {};
        this.lastScriptTaskResultPayload = resultPayload;
        const setText = (id, value) => {
            const el = document.getElementById(id);
            if (el) {
                el.textContent = value ?? '-';
            }
        };
        setText('script-task-result-task', resultPayload?.task_id || result.task_id || '-');
        setText('script-task-result-status', resultPayload?.job_status || '-');
        setText('script-task-result-exit-code', result.exit_code ?? '-');
        setText('script-task-output-dir', result.output_dir || '-');
        setText('script-task-artifact-count', Array.isArray(result.artifacts) ? result.artifacts.length : '-');
        if (this.scriptTaskResultPreview) {
            this.scriptTaskResultPreview.textContent = JSON.stringify({
                command: result.command || null,
                stdout_tail: result.stdout_tail || '',
                stderr_tail: result.stderr_tail || '',
                artifacts: result.artifacts || [],
                manifest_paths: result.manifest_paths || [],
            }, null, 2);
        }
        if (this.scriptTaskResultPanel) {
            this.scriptTaskResultPanel.style.display = '';
        }
    }

    async handleScriptTaskSubmit(event) {
        event.preventDefault();
        if (this.scriptTaskRunning) {
            return;
        }
        this.clearError();
        if (this.scriptTaskResultPanel) {
            this.scriptTaskResultPanel.style.display = 'none';
        }
        try {
            const task = this.getScriptTaskById();
            const payload = this.collectScriptTaskPayload();
            this.lastScriptTaskPayload = payload;
            if (this.scriptTaskCopyCurlBtn) {
                this.scriptTaskCopyCurlBtn.disabled = false;
            }
            const created = await API.createScriptTaskJob(payload, task);
            this.activeScriptTaskJobId = created.job_id;
            this.setScriptTaskStatus(`任务已创建: ${created.job_id}`, 'queued');
            this.writeScriptTaskDetails(created);
            upsert_task_history_entry({
                job_id: created.job_id,
                kind: created.kind,
                module: 'script-tasks',
                task_id: created.task_id,
                job_status: 'queued',
            });
            this.renderTaskHistory();
            this.startScriptTaskPolling(created.job_id);
        } catch (error) {
            this.setScriptTaskStatus(API.formatError(error), 'error');
            this.showError(API.formatError(error));
            this.showTechError(error);
        }
    }

    async handleScriptTaskCancel() {
        if (!this.activeScriptTaskJobId) {
            this.setScriptTaskStatus('当前没有可取消任务', 'error');
            return;
        }
        try {
            const payload = await API.cancelScriptTaskJob(this.activeScriptTaskJobId);
            this.stopScriptTaskPolling();
            this.setScriptTaskStatus(`任务已取消: ${payload.job_id}`, 'success');
            this.writeScriptTaskDetails(payload);
            upsert_task_history_entry({
                job_id: payload.job_id,
                kind: payload.kind,
                module: 'script-tasks',
                task_id: payload.task_id,
                job_status: payload.job_status,
            });
            this.renderTaskHistory();
        } catch (error) {
            this.setScriptTaskStatus(API.formatError(error), 'error');
            this.showTechError(error);
        }
    }

    async handleCopyScriptTaskCurl() {
        const text = this.buildScriptTaskCurlCommand();
        if (!text) {
            this.setScriptTaskStatus('当前没有可复制 curl 命令', 'error');
            return;
        }
        try {
            if (!navigator?.clipboard?.writeText) {
                throw new Error('clipboard API unavailable');
            }
            await navigator.clipboard.writeText(text);
            this.setScriptTaskStatus('已复制脚本任务 curl 命令', 'success');
        } catch (_error) {
            this.setScriptTaskStatus('复制失败，请从详情区手动复制', 'error');
        }
    }

    async loadScriptTaskFromHistory(jobId) {
        if (!jobId) {
            return;
        }
        this.activeScriptTaskJobId = jobId;
        this.switchServicePanel('script-tasks');
        try {
            const statusPayload = await API.getScriptTaskJobStatus(jobId);
            this.writeScriptTaskDetails(statusPayload);
            this.setScriptTaskStatus(`已加载脚本任务: ${jobId} (${statusPayload.job_status})`, 'info');
            if (statusPayload.task_id) {
                this.selectScriptTask(statusPayload.task_id);
            }
            if (['succeeded', 'failed', 'cancelled'].includes(statusPayload.job_status)) {
                const resultPayload = await API.getScriptTaskJobResult(jobId);
                this.updateScriptTaskResult(resultPayload);
            }
        } catch (error) {
            this.setScriptTaskStatus(API.formatError(error), 'error');
            this.showTechError(error);
        }
    }

    readNumberInput(id, fieldName) {
        const raw = String(document.getElementById(id)?.value || '').trim();
        if (!raw) {
            throw new Error(`${fieldName} 不能为空`);
        }
        const value = Number(raw);
        if (!Number.isFinite(value)) {
            throw new Error(`${fieldName} 必须是数字`);
        }
        return value;
    }

    readOptionalNumberInput(id, fieldName) {
        const raw = String(document.getElementById(id)?.value || '').trim();
        if (!raw) {
            return null;
        }
        const value = Number(raw);
        if (!Number.isFinite(value)) {
            throw new Error(`${fieldName} 必须是数字`);
        }
        return value;
    }

    readOptionalIntInput(id, fieldName) {
        const value = this.readOptionalNumberInput(id, fieldName);
        if (value === null) {
            return null;
        }
        if (!Number.isInteger(value)) {
            throw new Error(`${fieldName} 必须是整数`);
        }
        return value;
    }

    setPointStatus(statusEl, text, level = 'info') {
        if (!statusEl) {
            return;
        }
        statusEl.textContent = text;
        statusEl.className = `scan-status ${level}`;
    }

    setButtonBusy(button, busy, busyText, idleText) {
        if (!button) {
            return;
        }
        button.disabled = busy;
        button.textContent = busy ? busyText : idleText;
    }

    buildPointCurlCommand(url, payload) {
        if (!payload) {
            return '';
        }
        const body = JSON.stringify(payload);
        const escapedBody = body.replace(/'/g, "'\"'\"'");
        return `curl -X POST "${url}" -H "Content-Type: application/json" --data-raw '${escapedBody}'`;
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
        const mode = this.getScanMode();

        const payload = {
            kind,
            params: {
                mode,
                max_retries: 0,
            },
        };

        if (mode === 'point') {
            payload.params.T_mev = Number(document.getElementById('scan-point-t-mev')?.value);
            if (kind === 'trho') {
                payload.params.rho_value = Number(document.getElementById('scan-point-rho-value')?.value);
            } else {
                payload.params.mu_mev = Number(document.getElementById('scan-point-mu-mev')?.value);
            }
        } else {
            payload.params.T_values = this.parseNumberList(document.getElementById('scan-t-values')?.value, 'T_values');
            if (kind === 'trho') {
                payload.params.rho_values = this.parseNumberList(document.getElementById('scan-rho-values')?.value, 'rho_values');
            } else {
                payload.params.mu_values = this.parseNumberList(document.getElementById('scan-mu-values')?.value, 'mu_values');
            }
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

    initializeTemplates() {
        if (!this.scanTemplateSelect) {
            return;
        }
        const templates = load_task_templates();
        this.scanTemplateSelect.innerHTML = '';

        const defaultOption = document.createElement('option');
        defaultOption.value = '';
        defaultOption.textContent = '请选择模板';
        this.scanTemplateSelect.appendChild(defaultOption);

        const sourceTemplates = Array.isArray(templates) && templates.length > 0 ? templates : get_default_templates();
        sourceTemplates.forEach((tpl) => {
            const option = document.createElement('option');
            option.value = tpl.id;
            option.textContent = tpl.name;
            option.dataset.payload = JSON.stringify(tpl.payload);
            this.scanTemplateSelect.appendChild(option);
        });
    }

    applyTemplateSelection() {
        if (!this.scanTemplateSelect || !this.scanTemplateSelect.value) {
            return;
        }
        const selected = this.scanTemplateSelect.options[this.scanTemplateSelect.selectedIndex];
        const payload = selected?.dataset?.payload ? JSON.parse(selected.dataset.payload) : null;
        if (!payload) {
            return;
        }

        const kindEl = document.getElementById('scan-kind');
        if (kindEl && payload.kind) {
            kindEl.value = payload.kind;
        }

        const modeEl = document.getElementById('scan-mode');
        if (modeEl && payload?.params?.mode) {
            modeEl.value = payload.params.mode;
        }
        this.syncScanFormRows();

        const params = payload.params || {};
        const setText = (id, value) => {
            const el = document.getElementById(id);
            if (el) {
                el.value = value;
            }
        };

        if (Array.isArray(params.T_values)) {
            setText('scan-t-values', params.T_values.join(','));
        }
        if (params.T_mev !== undefined) {
            setText('scan-point-t-mev', String(params.T_mev));
        }
        if (Array.isArray(params.mu_values)) {
            setText('scan-mu-values', params.mu_values.join(','));
        }
        if (params.mu_mev !== undefined) {
            setText('scan-point-mu-mev', String(params.mu_mev));
        }
        if (Array.isArray(params.rho_values)) {
            setText('scan-rho-values', params.rho_values.join(','));
        }
        if (params.rho_value !== undefined) {
            setText('scan-point-rho-value', String(params.rho_value));
        }

        if (params.xi !== undefined) {
            const radio = document.querySelector('input[name="xi-strategy"][value="xi"]');
            if (radio) {
                radio.checked = true;
            }
            setText('scan-xi', String(params.xi));
        } else if (Array.isArray(params.xi_values)) {
            const radio = document.querySelector('input[name="xi-strategy"][value="xi_values"]');
            if (radio) {
                radio.checked = true;
            }
            setText('scan-xi-values', params.xi_values.join(','));
        } else if (params.xi_grid) {
            const radio = document.querySelector('input[name="xi-strategy"][value="xi_grid"]');
            if (radio) {
                radio.checked = true;
            }
            setText('scan-xi-grid-start', params.xi_grid.start);
            setText('scan-xi-grid-stop', params.xi_grid.stop);
            setText('scan-xi-grid-step', params.xi_grid.step);
        }

        this.showSuccess(`已应用模板: ${selected.textContent}`);
    }

    renderTaskHistory() {
        if (!this.scanHistoryList) {
            return;
        }
        const history = load_task_history();
        this.scanHistoryList.innerHTML = '';
        if (!history || history.length === 0) {
            this.scanHistoryList.textContent = '暂无历史任务';
            return;
        }

        history.slice(0, 10).forEach((item) => {
            const row = document.createElement('button');
            row.type = 'button';
            row.className = 'scan-history-item';
            const label = item.task_id ? `${item.task_id}` : (item.kind || '-');
            row.textContent = `${item.job_id} | ${item.job_status || 'unknown'} | ${label}`;
            row.addEventListener('click', () => {
                if (item.module === 'script-tasks' || item.kind === 'script-task') {
                    this.loadScriptTaskFromHistory(item.job_id);
                } else {
                    this.loadTaskFromHistory(item.job_id);
                }
            });
            this.scanHistoryList.appendChild(row);
        });

        this.renderHistoryHeatmap(history);
    }

    renderHistoryHeatmap(historyItems) {
        if (!this.historyHeatmap) {
            return;
        }
        this.historyHeatmap.innerHTML = '';
        const statuses = Array.isArray(historyItems) ? historyItems.slice(0, 20) : [];
        for (let i = 0; i < 20; i += 1) {
            const item = statuses[i];
            const cell = document.createElement('div');
            cell.className = `scan-heat-cell ${item?.job_status || 'empty'}`;
            cell.title = item ? `${item.job_id} | ${item.job_status}` : 'empty';
            this.historyHeatmap.appendChild(cell);
        }
    }

    drawProgressChart(percent) {
        if (!this.progressChart) {
            return;
        }
        const canvas = this.progressChart;
        const ctx = canvas.getContext('2d');
        if (!ctx) {
            return;
        }
        const value = Number.isFinite(Number(percent)) ? Number(percent) : 0;
        const w = canvas.width;
        const h = canvas.height;

        ctx.clearRect(0, 0, w, h);
        ctx.fillStyle = '#f8fafc';
        ctx.fillRect(0, 0, w, h);

        ctx.strokeStyle = '#cbd5e1';
        ctx.strokeRect(8, h / 2 - 10, w - 16, 20);

        ctx.fillStyle = '#38bdf8';
        const barWidth = Math.max(0, Math.min(w - 16, (w - 16) * value / 100));
        ctx.fillRect(8, h / 2 - 10, barWidth, 20);

        ctx.fillStyle = '#0f172a';
        ctx.font = '12px monospace';
        ctx.fillText(`Progress ${value.toFixed(1)}%`, 8, h / 2 - 16);
    }

    async loadTaskFromHistory(jobId) {
        if (!jobId) {
            return;
        }
        this.activeJobId = jobId;
        if (this.scanJobIdInput) {
            this.scanJobIdInput.value = jobId;
        }
        this.updateTaskRoute(jobId);
        await this.fetchAndRenderTask(jobId);
    }

    async handleLoadTaskByJobId() {
        const jobId = String(this.scanJobIdInput?.value || '').trim();
        if (!jobId) {
            this.setScanStatus('请输入 job_id', 'error');
            return;
        }
        this.activeJobId = jobId;
        this.updateTaskRoute(jobId);
        await this.fetchAndRenderTask(jobId);
    }

    updateTaskRoute(jobId = this.activeJobId) {
        if (typeof window === 'undefined') {
            return;
        }
        const routeParams = new URLSearchParams();
        if (jobId) {
            routeParams.set('job_id', jobId);
        }
        const query = routeParams.toString();
        window.location.hash = query ? `task-center?${query}` : 'task-center';
    }

    navigateToPage(pageId) {
        if (typeof window === 'undefined') {
            return;
        }
        if (pageId === 'task-center') {
            this.updateTaskRoute(this.activeJobId);
            return;
        }
        window.location.hash = 'scattering';
    }

    applyNavigationRoute() {
        if (typeof window === 'undefined') {
            return;
        }
        const hash = String(window.location.hash || '').replace(/^#/, '').trim();
        if (hash.startsWith('task-center')) {
            this.setActivePage('task-center');
            this.applyTaskRoute();
        } else {
            this.setActivePage('scattering');
        }
    }

    setActivePage(pageId) {
        const scatteringPage = document.getElementById('page-scattering');
        const taskCenterPage = document.getElementById('page-task-center');
        if (scatteringPage) {
            scatteringPage.classList.toggle('active', pageId === 'scattering');
        }
        if (taskCenterPage) {
            taskCenterPage.classList.toggle('active', pageId === 'task-center');
        }
        if (this.navScatteringBtn) {
            this.navScatteringBtn.classList.toggle('active', pageId === 'scattering');
        }
        if (this.navTaskCenterBtn) {
            this.navTaskCenterBtn.classList.toggle('active', pageId === 'task-center');
        }
        if (pageId === 'scattering' && this.visualization) {
            setTimeout(() => this.visualization.onWindowResize(), 0);
        }
    }

    applyTaskRoute() {
        if (typeof window === 'undefined') {
            return;
        }
        const hash = String(window.location.hash || '').replace(/^#/, '').trim();
        if (!hash.startsWith('task-center')) {
            return;
        }

        const queryIndex = hash.indexOf('?');
        const routeQuery = queryIndex >= 0 ? hash.slice(queryIndex + 1) : '';
        const jobId = new URLSearchParams(routeQuery).get('job_id');
        if (!jobId) {
            return;
        }
        if (this.scanJobIdInput) {
            this.scanJobIdInput.value = jobId;
        }
        this.activeJobId = jobId;
        this.fetchAndRenderTask(jobId);
    }

    async fetchAndRenderTask(jobId) {
        try {
            const statusPayload = await API.getJobStatus(jobId);
            this.writeScanDetails(statusPayload);
            this.drawProgressChart(statusPayload?.progress?.percent ?? 0);
            this.setScanStatus(`已加载任务: ${jobId} (${statusPayload.job_status})`, 'info');
            upsert_task_history_entry({
                job_id: statusPayload.job_id,
                kind: statusPayload.kind,
                job_status: statusPayload.job_status,
            });
            this.renderTaskHistory();

            if (statusPayload.job_status === 'succeeded') {
                const resultPayload = await API.getJobResult(jobId);
                this.updateScanResult(resultPayload);
                this.updateResultPreview(resultPayload);
            }
        } catch (error) {
            this.setScanStatus(API.formatError(error), 'error');
            this.showTechError(error);
        }
    }

    updateResultPreview(resultPayload) {
        const result = resultPayload?.result || {};
        this.lastResultPayload = resultPayload;
        if (this.scanResultPreview) {
            this.scanResultPreview.textContent = JSON.stringify(
                {
                    stats: result.stats || {},
                    output_path: result.output_path || null,
                },
                null,
                2,
            );
        }
        if (this.downloadResultBtn) {
            this.downloadResultBtn.disabled = !result.output_path;
        }
        this.updateScanReproInfo();
    }

    handleDownloadResult() {
        const outputPath = this.lastResultPayload?.result?.output_path;
        if (!outputPath) {
            this.setScanStatus('当前结果无可导出 output_path', 'error');
            return;
        }
        this.setScanStatus(`导出路径: ${outputPath}`, 'success');
    }

    getScanPayloadJson() {
        if (!this.lastScanPayload) {
            return '';
        }
        return JSON.stringify(this.lastScanPayload, null, 2);
    }

    buildScanCurlCommand() {
        if (!this.lastScanPayload) {
            return '';
        }
        const body = JSON.stringify(this.lastScanPayload);
        const escapedBody = body.replace(/'/g, "'\"'\"'");
        return `curl -X POST "${API.buildJobCreateUrl()}" -H "Content-Type: application/json" --data-raw '${escapedBody}'`;
    }

    updateScanReproInfo() {
        const hasPayload = !!this.lastScanPayload;
        if (this.scanReproPanel) {
            this.scanReproPanel.style.display = hasPayload ? '' : 'none';
        }
        if (this.copyParamsBtn) {
            this.copyParamsBtn.disabled = !hasPayload;
        }
        if (this.copyCurlBtn) {
            this.copyCurlBtn.disabled = !hasPayload;
        }
        if (this.scanReproPreview) {
            this.scanReproPreview.textContent = hasPayload
                ? JSON.stringify(
                    {
                        params: this.lastScanPayload,
                        create_job_curl: this.buildScanCurlCommand(),
                    },
                    null,
                    2,
                )
                : '-';
        }
    }

    async copyTextToClipboard(text, successMessage) {
        if (!text) {
            this.setScanStatus('当前没有可复制内容', 'error');
            return;
        }
        try {
            if (!navigator?.clipboard?.writeText) {
                throw new Error('clipboard API unavailable');
            }
            await navigator.clipboard.writeText(text);
            this.setScanStatus(successMessage, 'success');
        } catch (_error) {
            this.setScanStatus('复制失败，请手动从复现信息区块复制', 'error');
        }
    }

    async handleCopyScanParams() {
        await this.copyTextToClipboard(this.getScanPayloadJson(), '已复制参数 JSON');
    }

    async handleCopyScanCurl() {
        await this.copyTextToClipboard(this.buildScanCurlCommand(), '已复制 curl 命令');
    }

    collectPnjlGapPayload() {
        const rhoTarget = this.readOptionalNumberInput('pnjl-gap-rho-target', 'rho_target');
        const params = {
            T_mev: this.readNumberInput('pnjl-gap-t-mev', 'T_mev'),
            xi: this.readNumberInput('pnjl-gap-xi', 'xi'),
            allow_seed_fallback: !!document.getElementById('pnjl-gap-allow-seed-fallback')?.checked,
        };
        const muMev = rhoTarget === null
            ? this.readNumberInput('pnjl-gap-mu-mev', 'mu_mev')
            : this.readOptionalNumberInput('pnjl-gap-mu-mev', 'mu_mev');
        const pNum = this.readOptionalIntInput('pnjl-gap-p-num', 'p_num');
        const tNum = this.readOptionalIntInput('pnjl-gap-t-num', 't_num');
        if (muMev !== null) {
            params.mu_mev = muMev;
        }
        if (rhoTarget !== null) {
            params.rho_target = rhoTarget;
        }
        if (pNum !== null) {
            params.p_num = pNum;
        }
        if (tNum !== null) {
            params.t_num = tNum;
        }
        return { params };
    }

    updatePnjlGapResult(payload) {
        const result = payload?.result || {};
        this.lastPnjlGapResultPayload = payload;
        document.getElementById('pnjl-gap-converged').textContent = String(result.converged ?? '-');
        document.getElementById('pnjl-gap-iterations').textContent = result.iterations ?? '-';
        document.getElementById('pnjl-gap-residual').textContent = result.residual_norm ?? '-';
        document.getElementById('pnjl-gap-pressure').textContent = result.pressure ?? '-';
        if (this.pnjlGapResultPreview) {
            this.pnjlGapResultPreview.textContent = JSON.stringify(payload, null, 2);
        }
        if (this.pnjlGapResultPanel) {
            this.pnjlGapResultPanel.style.display = '';
        }
        if (this.pnjlGapCopyParamsBtn) {
            this.pnjlGapCopyParamsBtn.disabled = !this.lastPnjlGapPayload;
        }
        if (this.pnjlGapCopyCurlBtn) {
            this.pnjlGapCopyCurlBtn.disabled = !this.lastPnjlGapPayload;
        }
    }

    async handlePnjlGapSubmit(event) {
        event.preventDefault();
        this.clearError();
        this.setButtonBusy(this.pnjlGapSubmitBtn, true, '运行中...', '运行单点');
        if (this.pnjlGapResultPanel) {
            this.pnjlGapResultPanel.style.display = 'none';
        }
        try {
            const payload = this.collectPnjlGapPayload();
            this.lastPnjlGapPayload = payload;
            this.setPointStatus(this.pnjlGapStatus, '运行中...', 'running');
            const resultPayload = await API.runPnjlGap(payload);
            this.updatePnjlGapResult(resultPayload);
            this.setPointStatus(this.pnjlGapStatus, '运行完成', 'success');
        } catch (error) {
            this.setPointStatus(this.pnjlGapStatus, API.formatError(error), 'error');
            this.showError(API.formatError(error));
        } finally {
            this.setButtonBusy(this.pnjlGapSubmitBtn, false, '运行中...', '运行单点');
        }
    }

    async handleCopyPnjlGapParams() {
        const text = this.lastPnjlGapPayload ? JSON.stringify(this.lastPnjlGapPayload, null, 2) : '';
        await this.copyTextToClipboard(text, '已复制 PNJL 参数 JSON');
    }

    async handleCopyPnjlGapCurl() {
        await this.copyTextToClipboard(
            this.buildPointCurlCommand(API.buildPnjlGapUrl(), this.lastPnjlGapPayload),
            '已复制 PNJL curl 命令',
        );
    }

    collectTransportPointPayload() {
        const params = {
            T_mev: this.readNumberInput('transport-point-t-mev', 'T_mev'),
            mu_mev: this.readNumberInput('transport-point-mu-mev', 'mu_mev'),
            xi: this.readNumberInput('transport-point-xi', 'xi'),
            tau: this.readNumberInput('transport-point-tau', 'tau'),
            compute_bulk: !!document.getElementById('transport-point-compute-bulk')?.checked,
        };
        const pNum = this.readOptionalIntInput('transport-point-p-num', 'p_num');
        const tNum = this.readOptionalIntInput('transport-point-t-num', 't_num');
        const pNodes = this.readOptionalIntInput('transport-point-p-nodes', 'transport.p_nodes');
        const pMax = this.readOptionalNumberInput('transport-point-p-max', 'transport.p_max');
        const cosNodes = this.readOptionalIntInput('transport-point-cos-nodes', 'transport.cos_nodes');
        if (pNum !== null) {
            params.p_num = pNum;
        }
        if (tNum !== null) {
            params.t_num = tNum;
        }
        const transport = {};
        if (pNodes !== null) {
            transport.p_nodes = pNodes;
        }
        if (pMax !== null) {
            transport.p_max = pMax;
        }
        if (cosNodes !== null) {
            transport.cos_nodes = cosNodes;
        }
        if (Object.keys(transport).length > 0) {
            params.transport = transport;
        }
        return { params };
    }

    updateTransportPointResult(payload) {
        const result = payload?.result || {};
        const transport = result.transport || {};
        this.lastTransportPointResultPayload = payload;
        document.getElementById('transport-point-converged').textContent = String(result.equilibrium?.converged ?? '-');
        document.getElementById('transport-point-eta').textContent = transport.eta ?? '-';
        document.getElementById('transport-point-sigma').textContent = transport.sigma ?? '-';
        document.getElementById('transport-point-zeta').textContent = transport.zeta ?? '-';
        if (this.transportPointResultPreview) {
            this.transportPointResultPreview.textContent = JSON.stringify(payload, null, 2);
        }
        if (this.transportPointResultPanel) {
            this.transportPointResultPanel.style.display = '';
        }
        if (this.transportPointCopyParamsBtn) {
            this.transportPointCopyParamsBtn.disabled = !this.lastTransportPointPayload;
        }
        if (this.transportPointCopyCurlBtn) {
            this.transportPointCopyCurlBtn.disabled = !this.lastTransportPointPayload;
        }
    }

    async handleTransportPointSubmit(event) {
        event.preventDefault();
        this.clearError();
        this.setButtonBusy(this.transportPointSubmitBtn, true, '运行中...', '运行单点');
        if (this.transportPointResultPanel) {
            this.transportPointResultPanel.style.display = 'none';
        }
        try {
            const payload = this.collectTransportPointPayload();
            this.lastTransportPointPayload = payload;
            this.setPointStatus(this.transportPointStatus, '运行中...', 'running');
            const resultPayload = await API.runTransportPoint(payload);
            this.updateTransportPointResult(resultPayload);
            this.setPointStatus(this.transportPointStatus, '运行完成', 'success');
        } catch (error) {
            this.setPointStatus(this.transportPointStatus, API.formatError(error), 'error');
            this.showError(API.formatError(error));
        } finally {
            this.setButtonBusy(this.transportPointSubmitBtn, false, '运行中...', '运行单点');
        }
    }

    async handleCopyTransportPointParams() {
        const text = this.lastTransportPointPayload ? JSON.stringify(this.lastTransportPointPayload, null, 2) : '';
        await this.copyTextToClipboard(text, '已复制 Transport 参数 JSON');
    }

    async handleCopyTransportPointCurl() {
        await this.copyTextToClipboard(
            this.buildPointCurlCommand(API.buildTransportPointUrl(), this.lastTransportPointPayload),
            '已复制 Transport curl 命令',
        );
    }

    async handleScanCancel() {
        if (!this.activeJobId) {
            this.setScanStatus('当前没有可取消任务', 'error');
            return;
        }
        try {
            const payload = await API.cancelJob(this.activeJobId);
            this.stopPolling();
            this.setScanStatus(`任务已取消: ${payload.job_id}`, 'success');
            this.writeScanDetails(payload);
            upsert_task_history_entry({
                job_id: payload.job_id,
                kind: payload.kind,
                job_status: payload.job_status,
            });
            this.renderTaskHistory();
        } catch (error) {
            this.setScanStatus(API.formatError(error), 'error');
            this.showTechError(error);
        }
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
        this.updateScanReproInfo();
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
        if (this.scanCancelBtn) {
            this.scanCancelBtn.disabled = true;
        }
    }

    startPolling(jobId) {
        this.stopPolling();
        this.scanRunning = true;
        if (this.scanSubmitBtn) {
            this.scanSubmitBtn.disabled = true;
            this.scanSubmitBtn.textContent = '轮询中...';
        }
        if (this.scanCancelBtn) {
            this.scanCancelBtn.disabled = false;
        }

        this.pollTimer = setInterval(async () => {
            try {
                const statusPayload = await API.getJobStatus(jobId);
                this.writeScanDetails(statusPayload);
                const jobStatus = statusPayload?.job_status;
                const percent = statusPayload?.progress?.percent;
                this.drawProgressChart(percent ?? 0);
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
                } else if (jobStatus === 'cancelled') {
                    this.stopPolling();
                    this.setScanStatus('任务已取消', 'info');
                }

                upsert_task_history_entry({
                    job_id: statusPayload.job_id,
                    kind: statusPayload.kind,
                    job_status: statusPayload.job_status,
                });
                this.renderTaskHistory();
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
        this.lastResultPayload = null;
        try {
            const payload = this.collectScanPayload();
            this.lastScanPayload = payload;
            this.updateScanReproInfo();
            const created = await API.createScanJob(payload);
            this.activeJobId = created.job_id;
            this.setScanStatus(`任务已创建: ${created.job_id}`, 'queued');
            this.writeScanDetails(created);
            this.updateTaskRoute(created.job_id);
            if (this.scanJobIdInput) {
                this.scanJobIdInput.value = created.job_id;
            }
            upsert_task_history_entry({
                job_id: created.job_id,
                kind: created.kind,
                job_status: 'queued',
            });
            this.renderTaskHistory();
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
        this.lastResultPayload = null;
        this.updateScanReproInfo();
        try {
            const created = await API.createScanJob(this.lastScanPayload);
            this.activeJobId = created.job_id;
            this.setScanStatus(`重试已创建任务: ${created.job_id}`, 'queued');
            this.writeScanDetails(created);
            this.updateTaskRoute(created.job_id);
            if (this.scanJobIdInput) {
                this.scanJobIdInput.value = created.job_id;
            }
            upsert_task_history_entry({
                job_id: created.job_id,
                kind: created.kind,
                job_status: 'queued',
            });
            this.renderTaskHistory();
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
