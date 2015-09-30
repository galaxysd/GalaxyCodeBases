package main

//Usage: vsina http://video.weibo.com/show?fid=1034:0e906d53094c5d231bf09028af8ba9b1
import (
	"bufio"
	"fmt"
	"github.com/astaxie/beego/httplib"
	"io"
	"io/ioutil"
	"net/url"
	"os"
	"regexp"
	"strings"
	"time"
)

func help() {
	fmt.Println("Usage: vsina <VideoUrl>")
}

func fetch(videoUrl string) {
	req := httplib.Get(videoUrl)
	resp, err := req.Response()
	if err != nil {
		fmt.Println("fetch video url error,", err)
		return
	}
	defer resp.Body.Close()
	respData, err := ioutil.ReadAll(resp.Body)
	if err != nil {
		fmt.Println("read video url content error,", err)
		return
	}
	pattern := `flashvars="list=.*"`
	regx := regexp.MustCompile(pattern)
	flashVars := regx.FindString(string(respData))
	if flashVars != "" {
		size := len(flashVars)
		m3u8Url, err := url.QueryUnescape(flashVars[16 : size-1])
		if err != nil {
			fmt.Println(err)
		} else {
			fetchMovie(m3u8Url)
		}
	} else {
		fmt.Println("m3u8 playlist not found")
	}
}

func fetchMovie(m3u8Url string) {
	req := httplib.Get(m3u8Url)
	resp, respErr := req.Response()
	if respErr != nil {
		fmt.Println("fetch m3u8 playlist error,", respErr)
		return
	}
	defer resp.Body.Close()
	respData, readErr := ioutil.ReadAll(resp.Body)
	if readErr != nil {
		fmt.Println("read m3u8 playlist content error,", readErr)
		return
	}
	videoFolder := fmt.Sprintf("video_%d", time.Now().Unix())
	mkdErr := os.Mkdir(videoFolder, 0775)
	if mkdErr != nil {
		fmt.Println("mkdir for m3u8 playlist failed,", mkdErr)
		return
	}
	m3u8Uri, _ := url.Parse(m3u8Url)

	sReader := strings.NewReader(string(respData))
	bReader := bufio.NewScanner(sReader)
	bReader.Split(bufio.ScanLines)
	for bReader.Scan() {
		line := bReader.Text()
		if !strings.HasPrefix(line, "#") {
			tsFileName := line
			tsLocalFileName := fmt.Sprintf("%s/%s", videoFolder, line)
			tsFileUrl := fmt.Sprintf("%s://%s/%s", m3u8Uri.Scheme, m3u8Uri.Host, tsFileName)
			downloadTS(tsFileUrl, tsLocalFileName)
		}
	}
	fmt.Println()
	fmt.Println("Result:", videoFolder)
}

func downloadTS(tsFileUrl string, tsLocalFileName string) {
	fmt.Println("downloading", tsFileUrl)
	req := httplib.Get(tsFileUrl)
	resp, respErr := req.Response()
	if respErr != nil {
		fmt.Println("download ts ", tsFileUrl, "failed,", respErr)
		return
	}
	defer resp.Body.Close()
	tsFp, openErr := os.OpenFile(tsLocalFileName, os.O_CREATE|os.O_WRONLY, 0775)
	if openErr != nil {
		fmt.Println("open local file", tsLocalFileName, "failed,", openErr)
		return
	}
	defer tsFp.Close()
	_, copyErr := io.Copy(tsFp, resp.Body)
	if copyErr != nil {
		fmt.Println("download ts", tsFileUrl, " failed,", copyErr)
	}
}

func main() {
	argv := os.Args
	argc := len(argv)
	if argc != 2 {
		help()
		return
	}

	videoUrl := argv[1]
	fetch(videoUrl)
}
