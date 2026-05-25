**下载 YouTube Music 音乐**

```bash
yt-dlp "https://music.youtube.com/playlist?list=PLwqG5QrBfRXlq3fMjMQD3ga-Z7wKPw7FD" --yes-playlist -x --audio-format mp3 --audio-quality 0 --embed-metadata --embed-thumbnail --convert-thumbnails jpg -o "~/music/%(title)s - %(artist,uploader|Unknown)s.%(ext)s" --download-archive "~/music/downloaded.txt" --no-overwrites --downloader aria2c --downloader-args "aria2c:-x 16 -s 16 -k 1M"
```

- `--yes-playlist` 确认下载链接中的整个播放列表，而不仅仅是单首歌曲
- `-x` 提取音频（Extract audio），下载后自动剥离视频流，只保留声音
- `--audio-format mp3` 指定转换后的音频格式为 MP3
- `--audio-quality 0` 设置音频质量等级，`0` 代表最高质量（对于 MP3 是最好的可变码率）
- `--embed-metadata` 将歌曲标题、艺术家、专辑等信息写入文件标签（ID3 Tag）
- `--embed-thumbnail` 将视频封面图片嵌入到音频文件的封面中
- `--convert-thumbnails jpg` 将下载的封面图统一转换为 JPG 格式，以确保更好的兼容性
- `-o` 代表 output，定义输出路径和文件名模板（这里包含：标题 - 艺术家.扩展名）
- `--download-archive` 指定一个记录文件，记录已下载过的 ID，下次运行时自动跳过，避免重复下载
- `--no-overwrites` 如果本地已存在同名文件，则不再覆盖下载
- `--downloader aria2c` 调用多线程下载器 `aria2c` 来代替内置下载引擎，极大提升下载速度

**下载 Prohub 视频**

```bash
yt-dlp -f "bv*[ext=mp4]+ba[ext=m4a]/b[ext=mp4] / bv*+ba/b" --merge-output-format mp4 --no-warnings --retries infinite --fragment-retries infinite --ignore-errors --cookies-from-browser chrome -o "%(title)s.%(ext)s" -a url.txt
```

- `-f "bv*+ba/b"` 格式选择优先下载“最佳视频（bv*）”和“最佳音频（ba）”并合并，若不可用，则下载包含音视频的“最佳单文件（b）”
- `--merge-output-format mp4` 合并输出格式,，如果视频和音频流是分开下载的，最后统一封装（Mux）为 MP4 格式的文件
- `-a url.txt` 批量读取，`-a` 代表 batch-file，即从指定的文本文件（url.txt）中逐行读取所有链接并进行批量下载

**rclone**

查看 R2 存储桶列表：

```bash
rclone lsd r2:
```

新建存储桶：

```bash
rclone mkdir r2:sync
```

同步文件夹到 R2 存储桶：

```bash
rclone sync /home/jiaoyuan/desktop/sync r2:sync --transfers 16 --checkers 32 --s3-chunk-size 64M --s3-upload-concurrency 8 --fast-list -P
```

将 R2 存储桶的文件同步到本地：

```bash
rclone sync r2:sync /home/jiaoyuan/desktop/sync --transfers 16 --checkers 32 --s3-chunk-size 64M --fast-list -P
```

- `--transfers 16`，同时上传的文件数量（并发数）
- `--checkers 32`，扫描和对比本地与云端文件差异的并发数
- `--s3-chunk-size 64M`，大文件分块上传时每一块的大小（越大通常速度越快，但更占内存）
- `--s3-upload-concurrency 8`，单个文件内部同时上传的分块数量
- `--fast-list`，优化文件列表获取方式，减少 API 调用并大幅加快扫描速度
- `-P`，显示实时传输进度、速度、已用时间和预计剩余时间
- `sync`，同步操作，使云端目录与本地完全一致（注意：云端多出的文件会被删除）

可以将同步的命令设置为 Bash 别名便于终端直接调用：

```bash
r2up() {
   rclone sync /home/jiaoyuan/desktop/sync r2:sync --transfers 16 --checkers 32 --s3-chunk-size 64M --s3-upload-concurrency 8 --fast-list -P "$@"  
}
  
r2down() {
   rclone sync r2:sync /home/jiaoyuan/desktop/sync --transfers 16 --checkers 32 --s3-chunk-size 64M --fast-list -P "$@"  
}
```

挂载 R2 存储桶：

```bash
mkdir ~/r2 # 创建挂载点
rclone mount r2:sync ~/sync --vfs-cache-mode full &
```

取消挂载：

```bash
fusermount -u ~/r2
```

查看存储桶目录列表：

```bash
rclone ls r2:sync
```

删除存储桶：

```bash
rclone rmdir r2:sync # 删除空存储桶
rclone purge r2:sync # 删除非空存储桶
rclone delete r2:sync # 清空存储桶
```

- `--dry-run` 参数可以查看要被删除的目录列表

**查看 ArchLinux 滚动更新的次数**

```bash
echo "Since $(head -n1 /var/log/pacman.log | cut -d ' ' -f 1,2): $(grep -c 'full system upgrade' /var/log/pacman.log) full system upgrades"
```

**查看 ArchLinux 错误日志**

```bash
journalctl -b -p 0..4
```

- `-b` 代表 boot，只显示从本次开机以来生成的日志
- `-p` 代表 priority（优先级），用来过滤日志的严重程度，`0..4`: 指定优先级的范围，数字越小越严重，`0-3`: 各种级别的错误（紧急、警报、关键、普通错误），`4`: 

**传输文件到曙光超算**

```bash
/data/projects_data/rayfile-c -a harbin02.hpccube.com -P 65012 -u jiaoyuan -w b1ed11352f4690ce95-0d1c-43f5-888c-d88d0372270e -tm -no-meta -symbolic-links follow -retry 10 -retrytimeout 30 -o upload -d /02.crossover -s /home/  
jiaoyuan/desktop/all.hap2.chr.snp.gz
```