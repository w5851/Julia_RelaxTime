/**
 * Three.js 可视化模块
 * 负责3D场景渲染和动量椭球可视化
 */

import * as THREE from 'https://cdn.jsdelivr.net/npm/three@0.160.0/build/three.module.js';
import { OrbitControls } from 'https://cdn.jsdelivr.net/npm/three@0.160.0/examples/jsm/controls/OrbitControls.js';

export class Visualization {
    constructor(containerId) {
        this.container = document.getElementById(containerId);
        this.scene = null;
        this.camera = null;
        this.renderer = null;
        this.controls = null;
        
        // 场景对象
        this.arrows = {
            p1: null,
            p2: null,
            p3: null,
            p4: null
        };
        this.ellipsoid = null;
        this.axes = null;
        
        this.initialize();
    }

    /**
     * 初始化Three.js场景
     */
    initialize() {
        // 创建场景
        this.scene = new THREE.Scene();
        this.scene.background = new THREE.Color(0x1a1a2e);
        
        // 创建相机
        const aspect = this.container.clientWidth / this.container.clientHeight;
        this.camera = new THREE.PerspectiveCamera(60, aspect, 0.1, 1000);
        this.camera.position.set(10, 10, 10);
        this.camera.lookAt(0, 0, 0);
        
        // 创建渲染器
        this.renderer = new THREE.WebGLRenderer({ antialias: true });
        this.renderer.setSize(this.container.clientWidth, this.container.clientHeight);
        this.renderer.setPixelRatio(window.devicePixelRatio);
        this.container.appendChild(this.renderer.domElement);
        
        // 添加轨道控制器
        this.controls = new OrbitControls(this.camera, this.renderer.domElement);
        this.controls.enableDamping = true;
        this.controls.dampingFactor = 0.05;
        this.controls.minDistance = 5;
        this.controls.maxDistance = 50;
        
        // 添加光源
        this.addLights();
        
        // 添加坐标轴
        this.addCoordinateAxes();
        
        // 添加网格
        this.addGrid();
        
        // 监听窗口大小变化
        window.addEventListener('resize', () => this.onWindowResize());
        
        // 开始渲染循环
        this.animate();
    }

    /**
     * 添加光源
     */
    addLights() {
        // 环境光
        const ambientLight = new THREE.AmbientLight(0xffffff, 0.6);
        this.scene.add(ambientLight);
        
        // 主方向光
        const mainLight = new THREE.DirectionalLight(0xffffff, 0.8);
        mainLight.position.set(10, 10, 10);
        this.scene.add(mainLight);
        
        // 补光
        const fillLight = new THREE.DirectionalLight(0xffffff, 0.3);
        fillLight.position.set(-10, -10, -10);
        this.scene.add(fillLight);
    }

    /**
     * 添加坐标轴
     */
    addCoordinateAxes() {
        const axesHelper = new THREE.AxesHelper(5);
        this.axes = axesHelper;
        this.scene.add(axesHelper);
        
        // 添加轴标签
        this.addAxisLabels();
    }

    /**
     * 添加轴标签
     */
    addAxisLabels() {
        const createTextSprite = (text, color) => {
            const canvas = document.createElement('canvas');
            canvas.width = 256;
            canvas.height = 256;
            const context = canvas.getContext('2d');
            
            context.fillStyle = color;
            context.font = 'Bold 100px Arial';
            context.textAlign = 'center';
            context.textBaseline = 'middle';
            context.fillText(text, 128, 128);
            
            const texture = new THREE.CanvasTexture(canvas);
            const spriteMaterial = new THREE.SpriteMaterial({ map: texture });
            const sprite = new THREE.Sprite(spriteMaterial);
            sprite.scale.set(0.8, 0.8, 1);
            
            return sprite;
        };
        
        const xLabel = createTextSprite('X', '#ff0000');
        xLabel.position.set(5.5, 0, 0);
        this.scene.add(xLabel);
        
        const yLabel = createTextSprite('Y', '#00ff00');
        yLabel.position.set(0, 5.5, 0);
        this.scene.add(yLabel);
        
        const zLabel = createTextSprite('Z', '#0000ff');
        zLabel.position.set(0, 0, 5.5);
        this.scene.add(zLabel);
    }

    /**
     * 添加网格
     */
    addGrid() {
        const gridHelper = new THREE.GridHelper(20, 20, 0x444444, 0x222222);
        this.scene.add(gridHelper);
    }

    /**
     * 创建动量箭头
     * @param {Array<number>} vector - 动量向量
     * @param {number} color - 颜色
     * @returns {THREE.ArrowHelper} 箭头对象
     */
    createArrow(vector, color) {
        const origin = new THREE.Vector3(0, 0, 0);
        const direction = new THREE.Vector3(...vector).normalize();
        const length = Math.sqrt(vector[0]**2 + vector[1]**2 + vector[2]**2);
        const headLength = Math.min(length * 0.2, 0.5);
        const headWidth = headLength * 0.6;
        
        const arrow = new THREE.ArrowHelper(
            direction,
            origin,
            length,
            color,
            headLength,
            headWidth
        );
        
        return arrow;
    }

    /**
     * 创建椭球网格
     * @param {Object} ellipsoidParams - 椭球参数
     * @param {Array<number>} ellipsoidParams.center - 中心位置
     * @param {Array<Array<number>>} ellipsoidParams.axes_directions - 主轴方向
     * @param {Array<number>} ellipsoidParams.half_lengths - 半轴长
     * @returns {THREE.Mesh} 椭球网格
     */
    createEllipsoid(ellipsoidParams) {
        const { center, axes_directions, half_lengths } = ellipsoidParams;
        
        // 创建单位球体几何体
        const geometry = new THREE.SphereGeometry(1, 64, 32);
        
        // 创建半透明材质
        const material = new THREE.MeshPhongMaterial({
            color: 0x00c8ff,
            transparent: true,
            opacity: 0.3,
            side: THREE.DoubleSide,
            shininess: 100
        });
        
        const mesh = new THREE.Mesh(geometry, material);
        
        // 应用缩放（半轴长）
        mesh.scale.set(...half_lengths);
        
        // 应用旋转（主轴方向）
        // axes_directions是列向量，需要转置为行向量
        const rotationMatrix = new THREE.Matrix4();
        rotationMatrix.set(
            axes_directions[0][0], axes_directions[1][0], axes_directions[2][0], 0,
            axes_directions[0][1], axes_directions[1][1], axes_directions[2][1], 0,
            axes_directions[0][2], axes_directions[1][2], axes_directions[2][2], 0,
            0, 0, 0, 1
        );
        mesh.applyMatrix4(rotationMatrix);
        
        // 应用平移（中心位置）
        mesh.position.set(...center);
        
        // 添加线框
        const wireframe = new THREE.WireframeGeometry(geometry);
        const lineMaterial = new THREE.LineBasicMaterial({ 
            color: 0x00c8ff, 
            transparent: true, 
            opacity: 0.6 
        });
        const wireframeMesh = new THREE.LineSegments(wireframe, lineMaterial);
        wireframeMesh.scale.copy(mesh.scale);
        wireframeMesh.applyMatrix4(rotationMatrix);
        wireframeMesh.position.copy(mesh.position);
        
        // 创建组
        const group = new THREE.Group();
        group.add(mesh);
        group.add(wireframeMesh);
        
        return group;
    }

    /**
     * 更新可视化
     * @param {Object} data - 计算结果数据
     */
    update(data) {
        // 清除旧对象
        this.clearSceneObjects();
        
        const { momenta, ellipsoid } = data;
        
        // 创建动量箭头
        this.arrows.p1 = this.createArrow(momenta.p1, 0xff4444);
        this.arrows.p2 = this.createArrow(momenta.p2, 0x4444ff);
        this.arrows.p3 = this.createArrow(momenta.p3, 0x44ff44);
        this.arrows.p4 = this.createArrow(momenta.p4, 0xffaa00);
        
        this.scene.add(this.arrows.p1);
        this.scene.add(this.arrows.p2);
        this.scene.add(this.arrows.p3);
        this.scene.add(this.arrows.p4);
        
        // 创建椭球
        this.ellipsoid = this.createEllipsoid(ellipsoid);
        this.scene.add(this.ellipsoid);
        
        // 调整相机视角
        this.fitCameraToScene();
    }

    /**
     * 清除场景中的动态对象
     */
    clearSceneObjects() {
        // 移除箭头
        Object.values(this.arrows).forEach(arrow => {
            if (arrow) {
                this.scene.remove(arrow);
            }
        });
        
        // 移除椭球
        if (this.ellipsoid) {
            this.scene.remove(this.ellipsoid);
        }
        
        // 重置引用
        this.arrows = { p1: null, p2: null, p3: null, p4: null };
        this.ellipsoid = null;
    }

    /**
     * 自动调整相机以适应场景
     */
    fitCameraToScene() {
        // 计算场景边界
        const box = new THREE.Box3();
        
        if (this.ellipsoid) {
            box.expandByObject(this.ellipsoid);
        }
        
        Object.values(this.arrows).forEach(arrow => {
            if (arrow) {
                box.expandByObject(arrow);
            }
        });
        
        const center = box.getCenter(new THREE.Vector3());
        const size = box.getSize(new THREE.Vector3());
        const maxDim = Math.max(size.x, size.y, size.z);
        const fov = this.camera.fov * (Math.PI / 180);
        let cameraZ = Math.abs(maxDim / Math.tan(fov / 2)) * 1.5;
        
        this.camera.position.set(center.x + cameraZ, center.y + cameraZ, center.z + cameraZ);
        this.camera.lookAt(center);
        this.controls.target.copy(center);
        this.controls.update();
    }

    /**
     * 窗口大小变化处理
     */
    onWindowResize() {
        const width = this.container.clientWidth;
        const height = this.container.clientHeight;
        
        this.camera.aspect = width / height;
        this.camera.updateProjectionMatrix();
        
        this.renderer.setSize(width, height);
    }

    /**
     * 动画循环
     */
    animate() {
        requestAnimationFrame(() => this.animate());
        
        this.controls.update();
        this.renderer.render(this.scene, this.camera);
    }

    /**
     * 清理资源
     */
    dispose() {
        this.clearSceneObjects();
        this.renderer.dispose();
        this.controls.dispose();
        window.removeEventListener('resize', this.onWindowResize);
    }
}
